cd /doyale3/twarren/grebe2/r95pk42/wnon3
cpscript('/doyale/twarren/r95pk42/wnon3',42,01)
cpscript('/doyale/twarren/r95pk42/wnon3',42,02)

    cd /doyale/twarren/r95pk42/wnon3
    ls *42_01*_0*cbin>batch01a
    ls *42_01*_1*cbin>batch01b
    
    ls *42_02*_0*cbin>batch02a
    ls *42_02*_1*cbin>batch02b
    
    cleandir4tw('batch01a',1e4,500,6,10,'obs0')
    mk_rmdata('batch01a.dcrd',1)
    !csh rmdata ./
    cleandir4tw('batch01b',1e4,500,6,10,'obs0')
    mk_rmdata('batch01b.dcrd',1)
    !csh rmdata ./
    
   cleandir4tw('batch02a',1e4,500,6,10,'obs0')
    mk_rmdata('batch02a.dcrd',1)
    !csh rmdata ./
    
    cleandir4tw('batch02b',1e4,500,6,10,'obs0')
    mk_rmdata('batch02b.dcrd',1)
    !csh rmdata ./
   