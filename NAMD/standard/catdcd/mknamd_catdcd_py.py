#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("prefix", help="Prefix of the DCD")
parser.add_argument("--longname", action="store_true", help="It is a long name if containing pfXXps")
parser.add_argument("--ifreq", help="Frequency of input trajectory")
parser.add_argument("--ofreq", help="Frequency of output trajectory")
parser.add_argument("-b", "--begin", help="Begin index")
parser.add_argument("-e", "--end", help="End index")

args = parser.parse_args()

if args.longname:
    
