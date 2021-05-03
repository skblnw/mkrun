#!/bin/sed -rf
# Note: uses GNU sed extensions; might not work with other implementations

s/(^[ \t]*)\<fep\>/\1alchType           fep\n\1alch/I
s/(^[ \t]*)\<lambda\>/\1alchLambda/I
s/(^[ \t]*)\<lambda2\>/\1alchLambda2/I
s/(^[ \t]*)\<decouple\>/\1alchDecouple/I
/^[ \t]*\<fep[a-z]+\>/I { s/fep/alch/I }
