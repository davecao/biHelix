#### Download 4ZOW

    curl -s http://www.rcsb.org/pdb/files/4ZOW.pdb -o 4zow.pdb

#### SSE by stride

    mkdir -p stride
    ~/apps/stride/stride 4zow.pdb > stride/4zow.stride

#### Extract segments 

    mkdir -p segments;cd segments
    python ../extract_helix.py --pdbfile 4zow.pdb --pdbid 4ZOW --secstr stride/4zow.stride --ss-program=Stride

#### Calculate kinks
-k: calculate kinks only

    mkdir -p results;cd results

    for f in `ls segments/*.pdb | sed 's/.pdb//g'`;do fn=`basename $f`;../biHelix/bin/biHelix -i $f.pdb -q -k -o results/$fn._bihelix.pdb; done

