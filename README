# Installation inscructions:


first set your own $PIPELINE_HOME environment variable as the base of this package i.e. where the package has been checked out on your local system

1. Recompile bwa in $PIPELINE_HOME/soft/src/bwa* 

2. Recompile samtools in $PIPELINE_HOME/soft/src/samtools*
```shell
	cd $PIPELINE_HOME/soft/src/htslib
	./configure --prefix=$PIPELINE_HOME/soft/packages/htslib
	make
	make install
	cd $PIPELINE_HOME/soft/src/samtools
	./configure --prefix=$PIPELINE_HOME/soft/packages/samtools
	make
	make install
```
	
3. Recompile tabix in $PIPELINE_HOME/soft/src/tabix

4. Recompile verifyBAMId in $PIPELINE_HOME/soft/packages/verifyBamID/
```shell
	cd $PIPELINE_HOME/soft/packages/verifyBamID
	make cloneLib # do this if ../libStatGen does not already exist
	make
```

5. Recompile queueInterpreter in $PIPELINE_HOME/soft/src/queueInterpreter

6. Reference files not included in this repo. You must download them separately and manually create link to the reference files directory ($PIPELINE_HOME/data/reference)
    NOTE: all .bed, .vcf and related files used by pipeline scripts use GRCh37 coordinates. If you wish to use a different reference, be sure to modify all necessary values in data/templates/init_pipeline.sh, in addition to any user-defined config files and command arguments
	

# How to Launch pipeline:

run the main launcher script (launch_pipeline.MIPs.fastqs.scatter.sh) with " -h " or " --help " for a full list of mandatory and optional command-line arguments


to use individual functions from the function library (i.e. for testing, debugging etc)

```
export PIPELINE_HOME=< your pipeline root directory >
. $PIPELINE_HOME/data/templates/init_pipeline.sh $PIPELINE_HOME
. $PIPELINE_HOME/data/templates/function.library 
sample_init BA000VIL0001.test ./
[some function] [function arguments]
```