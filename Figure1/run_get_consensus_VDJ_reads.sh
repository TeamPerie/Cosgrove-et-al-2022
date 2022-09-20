#!/bin/bash
#  run_get_consensus_VDJ.sh
#
#  Created by A Lyne on 20/08/2018.
#

WD="$AGREP_DIR"
cd "$WD"

rm consensus_*
rm blast*
rm fasta*
rm VDJ_*

#R script to get consensus barcodes within CB-UMI and then within each cell
R_FILE="$SCRIPT_DIR/consensus_per_UMI_and_cell_reads.R"
$R_EXEC $R_FILE $WD $NMATCH_V $NMATCH_J $SAMPLE $MIN_READ $MIN_PROP $MIN_UMI $MIN_PROP_1UMI


#add flags to say if extracted VDJ barcodes conform to expected structure
R_FILE="$SCRIPT_DIR/splitVDJ.R"
$R_EXEC $R_FILE $WD $NMATCH_V $NMATCH_J $SAMPLE $MIN_READ $MIN_PROP $MIN_UMI $MIN_PROP_1UMI


