#!/bin/bash

file=${1:-MIP_QC.mip-by-mip}

awk -v c=0 '
BEGIN{FS=OFS="\t"; SUBSEP="@"}
{
    if (NR==1) next;
    if (NR==2) {
        for (i=1; i<=NF; i++) {
            if ($i=="onMED") callabilityLimit=i;
            if ($i=="mip") mipCol=i
        };
        for (i=mipCol-6; i<mipCol; i++) {
            COUNTS[$i]=1;
            countsCols[i]=$i
        };
        next
    }; 
    mip=$mipCol;
    mips[mip]=1;
    for (i=1; i<callabilityLimit; i++) {
        callability_data[$i,mip]++;
        STATES[$i]=1
    }
    for (i in countsCols) {
        countType=countsCols[i]
        count_data[countType,mip]=$i
    };
}
END{
    header="mip";
    asorti(STATES, sortedSTATES)
    asorti(COUNTS, sortedCOUNTS)
    asorti(mips, sortedMips)
    for (i in sortedSTATES) {
        state=sortedSTATES[i]
        if (state=="") state="NO_DATA";
        header=header"\t"state
    };
    for (i in sortedCOUNTS) {
        countType=sortedCOUNTS[i]
        header=header"\t"countType
    };
    print header;
    for (i in sortedMips) {
        mip=sortedMips[i];
        s="";
        c="";
        for (i in sortedSTATES) {
            state=sortedSTATES[i]
            if ((state,mip) in callability_data) {
                value=callability_data[state,mip]
            } else {
            value=0
        };
        s=s"\t"value
        };
        for (i in sortedCOUNTS) {
            countType=sortedCOUNTS[i]
            count=count_data[countType,mip]
            c=c"\t"count
        };
        printf "%s%s%s\n", mip,s,c
    };
}
' $file
