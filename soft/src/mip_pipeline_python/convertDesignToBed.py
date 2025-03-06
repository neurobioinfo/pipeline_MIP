#!/usr/bin/env python

import sys
import os
import Design

designFile, designSource = sys.argv[1:]
outputDir = os.path.dirname(designFile) + '/' if os.path.dirname(designFile) else ""
bedFile = outputDir + '.'.join(os.path.basename(designFile).split('.')[0:-1]) + '.bed'
design = Design.Design(designFile, designSource)
design.parse()
design.generateBedFiles(bedFile)
