if [ -z "${1}" ]
then
    export MAX_SEG=100000000000
else
    export MAX_SEG=${1}
fi

find traj_segs -name "seg.last_frame.aligned_to_ref_pocket.pdb" | awk -v MAX_SEG=${MAX_SEG} '{if (substr($0,11,6) + 0 <= MAX_SEG) {print "cat " $1 "; echo ENDMDL"}}' | parallel > vis/all_subpex.pdb 

