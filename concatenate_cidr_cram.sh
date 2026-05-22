#!/bin/bash
set -euxo pipefail

pref=/scratch/ucgd/lustre-core/UCGD_Research/quinlan_NIH/NIH_CIDR_CEPH/Quinlan_Released_Data/CRAM

module load samtools

samtools merge -O CRAM \
               -o ${pref}/300100.cram \
               -@ 8 \
                ${pref}/639719-3553507751.cram ${pref}/639544-3553507546.cram

samtools merge -O CRAM \
               -o ${pref}/80010.cram \
               -@ 8 \
                ${pref}/639640-3553499421.cram ${pref}/639484-3553500728.cram

samtools merge -O CRAM \
               -o ${pref}/NA12878.cram \
               -@ 8 \
                ${pref}/NA12878-3544045660.cram ${pref}/NA12878-3544045670.cram