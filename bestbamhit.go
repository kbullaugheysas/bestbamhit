package main

import (
	"compress/gzip"
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"log"
	"math/rand"
	"os"
	"strings"
	"time"
)

type Args struct {
	MinScore     int
	MaxDist      int
	Limit        int
	Penalty      float64
	Labels       string
	LogFilename  string
	KeepFilename string
}

type Hit struct {
	Index  int
	Record BamRecord
}

var args = Args{}
var logger *log.Logger

func init() {
	log.SetFlags(0)
	flag.IntVar(&args.MinScore, "min-score", 60, "min alignment score for an alignment")
	flag.IntVar(&args.MaxDist, "max-dist", 5, "max edit distance for an alignment")
	flag.IntVar(&args.Limit, "limit", 0, "limit the number of sample reads considered (0 = no limit)")
	flag.Float64Var(&args.Penalty, "edit-penalty", 2.0, "multiple for how to penalize edit distance")
	flag.StringVar(&args.Labels, "labels", "", "comma-separated list of labels for the BAMs (required)")
	flag.StringVar(&args.LogFilename, "log", "", "write parameters and stats to a log file")
	flag.StringVar(&args.KeepFilename, "keep", "", "file where to write the names of reads matching the first bam file")

	flag.Usage = func() {
		log.Println("usage: bestbamhit [options] a.bam b.bam ...")
		flag.PrintDefaults()
	}
}

func benchmark(start time.Time, label string) {
	elapsed := time.Since(start)
	logger.Printf("%s took %s", label, elapsed)
}

func OpenLogger() {
	if args.LogFilename == "" {
		logger = log.New(os.Stderr, "", 0)
	} else {
		logfile, err := os.Create(args.LogFilename)
		if err != nil {
			log.Fatal(err)
		}
		logger = log.New(logfile, "", 0)
	}
}

func LogArguments() {
	logger.Println("command:", strings.Join(os.Args, " "))
	blob, err := json.MarshalIndent(args, "", "    ")
	if err != nil {
		logger.Fatal("failed to marshal arguments")
	}
	logger.Println(string(blob))
}

func main() {
	flag.Parse()
	bams := flag.Args()
	startedAt := time.Now()

	OpenLogger()

	if len(bams) == 0 {
		logger.Println("must specify at least one BAM file")
		os.Exit(1)
	}

	if args.Labels == "" {
		logger.Println("must specify -labels lab1,lab2")
		os.Exit(1)
	}

	labels := strings.Split(args.Labels, ",")

	LogArguments()

	var err error
	var keep_writer io.Writer
	if args.KeepFilename != "" {
		keepfile, err := os.Create(args.KeepFilename)
		if err != nil {
			log.Fatal(err)
		}
		if strings.HasSuffix(args.KeepFilename, ".gz") {
			keepz := gzip.NewWriter(keepfile)
			defer func() {
				if err := keepz.Close(); err != nil {
					log.Fatal(err)
				}
			}()
			keep_writer = keepz
		} else {
			defer keepfile.Close()
			keep_writer = keepfile
		}
	}

	scanners := make([]*BamScanner, len(bams))

	for c := 0; c < len(bams); c++ {
		scanners[c], err = OpenBam(bams[c])
		if err != nil {
			logger.Fatal(err)
		}
	}

	total_mappings := 0
	too_diverged := 0
	tooLow := 0
	found := 0
	ercc := 0
	multi := 0
	totalMatchLen := 0
	accepted := 0
	counts := make([]int, len(bams))

	err = func() error {
		defer benchmark(startedAt, "processing")

		for {

			if found > 0 && found%100000 == 0 {
				logger.Printf("found %d reads so far\n", found)
			}
			if args.Limit > 0 && args.Limit == found {
				return nil
			}

			var read string

			// Find the indexes of the scanners that present the lowest-ordered read.
			all_closed := true
			for _, s := range scanners {
				record, err := s.Record()
				if err != nil {
					return err
				}
				if s.Closed {
					continue
				}
				all_closed = false
				if read == "" {
					read = record[0]
				} else {
					if strnum_cmp(record[0], read) < 0 {
						read = record[0]
					}
				}
			}
			if all_closed {
				return nil
			}
			if read == "" {
				return fmt.Errorf("Failed to find read")
			}

			// Get all the records from all the bams that are for this read.
			var hits []Hit
			for i, s := range scanners {
				for {
					record, err := s.Record()
					if err != nil {
						return err
					}
					if record == nil || record[0] != read {
						break
					}
					hit := Hit{Index: i}
					err = hit.Record.Load(record)
					if err != nil {
						return err
					}
					hits = append(hits, hit)
					s.Ratchet()
				}
			}
			if len(hits) == 0 {
				return fmt.Errorf("No hits for %s", read)
			}
			total_mappings += len(hits)

			// Determine which of the hits has the best alignment
			var best_score int
			var which_best []int
			for j, hit := range hits {
				score := hit.Record.TagAS
				if len(which_best) == 0 {
					best_score = score
					which_best = append(which_best, j)
				} else {
					if score == best_score {
						which_best = append(which_best, j)
					} else if score > best_score {
						best_score = score
						which_best = []int{j}
					}
				}
			}
			var best int
			// If there are multiple best hits then we randomly select one.
			if len(which_best) > 1 {
				best = which_best[rand.Intn(len(which_best))]
				first_source := hits[which_best[0]].Index

				// Check if they're not all identical, in which case increment multi
				for _, j := range which_best {
					if hits[j].Index != first_source {
						multi++
						break
					}
				}
			} else {
				best = which_best[0]
			}
			// Check if it's ERCC
			if strings.Contains(hits[best].Record.Rname, "ERCC") {
				ercc++
			} else {
				bestHit := hits[best]
				mlen := bestHit.Record.MatchLength
				if err != nil {
					return err
				}
				if bestHit.Record.TagnM > args.MaxDist {
					too_diverged++
				} else if bestHit.Record.TagAS < args.MinScore {
					tooLow++
				} else {
					totalMatchLen += mlen
					accepted += 1
					counts[bestHit.Index] += 1
					if keep_writer != nil {
						fmt.Fprintf(keep_writer, "%s\t%s\n", hits[best].Record.Qname, labels[bestHit.Index])
					}
				}
			}
			found++
		}
	}()
	if err != nil {
		logger.Fatal(err)
	}

	avgMatchLen := float64(totalMatchLen) / float64(accepted)

	logger.Printf("total\t%d\n", total_mappings)
	logger.Printf("score too low\t%d\n", tooLow)
	logger.Printf("too diverged\t%d\n", too_diverged)
	logger.Printf("reads\t%d\n", found)
	logger.Printf("ercc\t%d\n", ercc)
	logger.Printf("multi\t%d\n", multi)
	logger.Printf("average match length\t%f\n", avgMatchLen)

	stats := []int{total_mappings, tooLow, too_diverged, found, ercc, multi}
	for c, count := range counts {
		logger.Printf("%s\t%d\n", labels[c], count)
		stats = append(stats, count)
	}

	statsStr := "stats"
	for _, s := range stats {
		statsStr += fmt.Sprintf("\t%d", s)
	}
	logger.Println(statsStr)

}
