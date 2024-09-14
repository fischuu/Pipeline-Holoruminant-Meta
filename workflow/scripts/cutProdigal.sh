awk -v size=$1 -v pre=$2 -v pad=7 '
     /^>/ { n++; if (n % size == 1) { close(fname); fname = sprintf("%s.%0" pad "d", pre, n) } }
     { print >> fname }
' $3

#find $2* | 's#.*/##' > $2.txt