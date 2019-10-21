#!/bin/bash -i
#ssh anton2.psc.edu

garden load msys/1.7.213c7/bin
prefix=md
mae2dms ${prefix}_converted.cms ${prefix}_out.dms
