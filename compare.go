package main

import (
	"log"
	"strconv"
	"unicode"
)

func digitToInt(s string) int {
	n, err := strconv.Atoi(s)
	if err != nil {
		log.Fatal("failed to parse digit ", err)
	}
	return n
}

// From: https://github.com/samtools/samtools/blob/develop/bam_sort.c#L13
func strnum_cmp(as, bs string) int {
	a := []rune(as)
	b := []rune(bs)
	i := 0
	j := 0
	for i < len(a) && j < len(b) {
		if unicode.IsDigit(a[i]) && unicode.IsDigit(b[j]) {
			for i < len(a) && a[i] == '0' {
				i++
			}
			for j < len(b) && b[j] == '0' {
				j++
			}
			for i < len(a) && j < len(b) && unicode.IsDigit(a[i]) && unicode.IsDigit(b[j]) && a[i] == b[j] {
				i++
				j++
			}
			// By this point we've forwarded across any leading zeros && any digits that match.
			// Next we get determine if they have the same number of digits
			// before the first non-diget. If so we use the numerical values of
			// the number formed by these digits to determine order.
			if i < len(a) && j < len(b) && unicode.IsDigit(a[i]) && unicode.IsDigit(b[j]) {
				k := 0
				for i+k < len(a) && unicode.IsDigit(a[i+k]) && j+k < len(b) && unicode.IsDigit(b[j+k]) {
					k += 1
				}
				if i+k < len(a) && unicode.IsDigit(a[i+k]) {
					return 1
				} else if j+k < len(b) && unicode.IsDigit(b[j+k]) {
					return -1
				} else {
					return digitToInt(string(a[i:(i+k)])) - digitToInt(string(b[j:(j+k)]))
				}
			} else if i < len(a) && unicode.IsDigit(a[i]) {
				return 1
			} else if j < len(b) && unicode.IsDigit(b[j]) {
				return -1
			} else if i != j {
				if i < j {
					return 1
				}
				return -1
			}
		} else {
			if a[i] != b[j] {
				return digitToInt(string(a[i])) - digitToInt(string(b[j]))
			}
			i++
			j++
		}
	}
	if len(a) > len(b) {
		return 1
	} else if len(a) < len(b) {
		return -1
	}
	return 0
}
