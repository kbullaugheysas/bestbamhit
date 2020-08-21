package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"os/exec"
	"regexp"
	"strconv"
	"strings"
	"sync"
)

type BamScanner struct {
	LineNumber int
	filename   string
	stdin      bool
	scanner    *bufio.Scanner
	wg         sync.WaitGroup
	prev       string
	record     []string
	Closed     bool
}

type BamRecord struct {
	Qname       string
	Flag        int
	Rname       string
	Pos         int
	Mapq        int
	Cigar       string
	Rnext       string
	Pnext       int
	Seq         string
	Qual        string
	TagAS       int
	TagHI       int
	TagnM       int
	MatchLength int
}

var cigarPattern *regexp.Regexp = regexp.MustCompile(`[0-9][0-9]*[A-Z]`)

func (r *BamRecord) Load(record []string) error {
	r.Qname = record[0]
	r.Rname = record[2]
	r.Cigar = record[5]
	r.Rnext = record[6]
	r.Seq = record[9]
	r.Qual = record[10]
	var err error
	if r.Flag, err = strconv.Atoi(record[1]); err != nil {
		return err
	}
	if r.Pos, err = strconv.Atoi(record[3]); err != nil {
		return err
	}
	if r.Mapq, err = strconv.Atoi(record[4]); err != nil {
		return err
	}
	if r.Pnext, err = strconv.Atoi(record[7]); err != nil {
		return err
	}
	for i := 11; i < len(record); i++ {
		val, err := strconv.Atoi(record[i][5:])
		if err != nil {
			return fmt.Errorf("failed to parse tag: %s", record[i])
		}
		switch record[i][0:5] {
		case "AS:i:":
			r.TagAS = val
		case "HI:i:":
			r.TagHI = val
		case "nM:i:":
			r.TagnM = val
		}
	}

	matches := cigarPattern.FindAllStringSubmatch(r.Cigar, -1)
	for _, match := range matches {
		code := match[0]
		if len(code) > 0 && code[len(code)-1] == 'M' {
			n, err := strconv.Atoi(code[0:(len(code) - 1)])
			if err != nil {
				return fmt.Errorf("failed to parse cigar fragment: %s", code)
			}
			r.MatchLength += n
		}
	}

	return nil
}

func OpenBam(bamfile string) (*BamScanner, error) {
	s := BamScanner{}
	s.filename = bamfile
	cmd := exec.Command("samtools", "view", bamfile)
	input, err := cmd.StdoutPipe()
	if err != nil {
		return nil, fmt.Errorf("failed creating pipe: %v", err)
	}
	if err := cmd.Start(); err != nil {
		return nil, fmt.Errorf("command failed to start: %v", err)
	}
	s.scanner = bufio.NewScanner(input)
	s.wg.Add(1)
	go func() {
		s.wg.Wait()

		if !s.stdin {
			if err := cmd.Wait(); err != nil {
				log.Fatal("wait failed: ", err)
			}
		}
	}()
	return &s, nil
}

// Fast forward to the next record with read name `read`
func (s *BamScanner) Find(read string) ([]string, error) {
	for {
		// The end of the file may have been reached previously.
		if s.Closed {
			return nil, nil
		}
		record, err := s.Record()
		if err != nil {
			return nil, err
		}
		// Or maybe the file is only now realized to be at the end.
		if s.Closed {
			return nil, nil
		}
		if record[0] == read {
			s.Ratchet()
			return record, nil
		}
		if strnum_cmp(record[0], read) < 0 {
			// Not far enough yet
			s.Ratchet()
		} else {
			// We didn't find the read before we reached one that is past what
			// we're looking for. We'll leave this one in the cache in case we
			// search for it next.
			return nil, nil
		}
	}
}

func (s *BamScanner) Record() ([]string, error) {
	if s.record != nil {
		return s.record, nil
	}
	s.Closed = !s.scanner.Scan()
	if err := s.scanner.Err(); err != nil {
		return nil, fmt.Errorf("scanner of %s errored: %v", s.filename, err)
	}
	if s.Closed {
		return nil, nil
	}
	line := strings.TrimSpace(s.scanner.Text())
	s.LineNumber++
	if len(line) == 0 {
		return nil, fmt.Errorf("empty BAM record")
	}
	s.record = strings.Split(line, "\t")
	if len(s.record) == 0 {
		return nil, fmt.Errorf("empty record at line %s", s.LineNumber)
	}
	read := s.record[0]
	if s.prev != "" {
		if strnum_cmp(s.prev, read) > 0 {
			return nil, fmt.Errorf("sorting order violated at line %d", s.LineNumber)
		}
	}
	s.prev = read
	return s.record, nil
}

func (s *BamScanner) Ratchet() {
	s.record = nil
}

func (s *BamScanner) Done() {
	s.wg.Done()
}

func ReadBamHeader(bamfile string) (string, error) {
	output, err := exec.Command("samtools", "view", "-H", bamfile).Output()
	if err != nil {
		return "", fmt.Errorf("failed to read header: %v", err)
	}
	return string(output), nil
}

type BamWriter struct {
	filename string
	wg       sync.WaitGroup
	fp       *os.File
}

func (w *BamWriter) Wait() {
	w.wg.Wait()
}
