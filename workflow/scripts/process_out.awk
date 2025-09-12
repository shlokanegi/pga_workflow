#!/usr/bin/awk -f

BEGIN {
    FS="\t"
}

{
    if ($5 != "*") {
        print $0
    }
}