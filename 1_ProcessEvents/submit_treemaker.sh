#!/bin/bash
# Script to simplify job submission to the grid
#
# Originally from Steven Gardiner <gardiner@fnal.gov>

##### Parse any command-line options passed to this script ######
# This section is based on https://stackoverflow.com/a/29754866/4081973
#
# Currently, the allowed options for invoking this script are
#
# -p, --prod:  Submit jobs with the VOMS Production role (rather than the
#              usual Analysis role). The job will be run on the grid
#              using the anniepro account instead of your normal user
#              account. Job submission will fail unless you have permission
#              to submit production jobs.
NUM_EXPECTED=5
NUM_EXPECTED_DEFAULTS=5

GRID_RESOURCES_DIR=""
COPY_FILES=""


# Check for a suitable version of the getopt program (should be installed
# on just about any flavor of Linux these days)
getopt --test > /dev/null
if [[ $? -ne 4 ]]; then
    echo "ERROR: enhanced getopt is not installed on this system."
    exit 1
fi

# Define the command-line options accepted by this script
OPTIONS=p:f::
LONGOPTIONS=prod:,file::

# Parse the command line used to run this script using getopt
PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTIONS --name "$0" -- "$@")

# Parse comand-line options 

if [[ $? -ne 0 ]]; then
    # getopt has complained about wrong arguments to stdout, so we can quit
    exit 2
fi

# Read getopt’s output this way to handle the quoting right
eval set -- "$PARSED"

# Process the options in order until we see --
while true; do
    case "$1" in
        -p|--prod)
            USE_PRODUCTION="y"
            shift
            ;;
        -f|--file)
            shift
	    COPY_FILES+="-f $1"
            echo "Copying files to node: $COPY_FILES"
	    shift
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Error encountered while parsing the command-line options: $1"
            exit 3
            ;;
    esac
done

if [ "$#" == "$NUM_EXPECTED_DEFAULTS" ]; then
    echo "Using default git tags: GENERATOR_GIT_CHECKOUT = $GENERATOR_GIT_CHECKOUT_DEFAULT, REWEIGHT_GIT_CHECKOUT = $REWEIGHT_GIT_CHECKOUT_DEFAULT, NUISANCE_GIT_CHECKOUT = $NUISANCE_GIT_CHECKOUT_DEFAULT and default MEC spline directory: RW_SPLINES_FILE = $RW_SPLINES_FILE_DEFAULT"
    NUM_JOBS=$1
    STEM=$2
    GPREP_FILE=$3
    CARD_FILE=$4
    OUTPUT_DIRECTORY=$5



else
  echo "Usage:"
  echo "Option 1: ./submit_nuisance.sh [--prod] NUM_JOBS JOB_NAME GPREP_FILE CARD_FILE OUTPUT_DIRECTORY RW_SPLINES_FILE GENERATOR_GIT_CHECKOUT REWEIGHT_GIT_CHECKOUT NUISANCE_GIT_CHECKOUT"
  echo "Option 2 (will use default values for RW_SPLINES_FILE, GENERATOR_GIT_CHECKOUT, REWEIGHT_GIT_CHECKOUT, AND NUISANCE_GIT_CHECKOUT - defaults are set at the top of this script): ./submit_nuisance.sh [--prod] NUM_JOBS JOB_NAME GPREP_FILE CARD_FILE OUTPUT_DIRECTORY"
  exit 1
fi

# Check if the USE_PRODUCTION variable is set
if [ ! -z ${USE_PRODUCTION+x} ]; then
  echo "Production job will run as uboonepro"
  PROD_OPTION="--role=Production"
else
  echo "Analysis job will run as $(whoami)"
  PROD_OPTION="--role=Analysis"
fi

source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups
setup fife_utils

jobsub_submit -G uboone --disk=20GB --memory=6000MB --expected-lifetime=96h -N $NUM_JOBS \
  $PROD_OPTION -f $GPREP_FILE -f $CARD_FILE \
  $COPY_FILES \
  -d OUTPUT $OUTPUT_DIRECTORY \
  --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE \
  file://$GRID_RESOURCES_DIR/treemaker_grid.sh $STEM \
  $(basename $GPREP_FILE) $(basename $CARD_FILE) 
