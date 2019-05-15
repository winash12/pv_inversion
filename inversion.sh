#!/bin/csh

# -------------------------------------------------------------------
# Set some variables and paths
# -------------------------------------------------------------------

# Jump to specified step
set step="help"
set param1= 
set param2=
set param3=

if ( $#argv == 1 ) then
  set step=$1
endif
if ( $#argv == 2 ) then
  set step=$1
  set param1=$2
endif
if ( $#argv == 3 ) then
  set step=$1
  set param1=$2
  set param2=$3
endif
if ( $#argv == 4 ) then
  set step=$1
  set param1=$2
  set param2=$3
  set param3=$4
endif

set file1=${param1}
set file2=${param2}

# Set base directory for programmes
#set bdir=${DYN_TOOLS}/inversion/
set bdir=/home/winash12/potentialVorticity

# Set name of the parameter file, coastline file and sample directory
set parafile="${PWD}/inversion.param"
set coastfile="${bdir}/prep/coastline.dat"
set sampledir="/net/dansgaard/atmosdyn/erainterim/cdf/2006/01/"

# ---------------------------------------------------------------------------
# Installation, help and sample
# ---------------------------------------------------------------------------

# Get parameter file, if not yet available
if ( ! -f ${parafile} ) then
       cp ${bdir}/inversion.param . 
endif

# Extract parameters from parameter file
set progex=${bdir}/inversion.perl


set    date=`${progex} ${parafile} date           | awk '{ print $2}'`


set nofiter=`${progex} ${parafile} n_of_iteration | awk '{ print $2}'`
set    save=`${progex} ${parafile} save_iteration | awk '{ print $2}'`
set    idir=`${progex} ${parafile} inp_dir        | awk '{ print $2}'`
set    rdir=`${progex} ${parafile} run_dir        | awk '{ print $2}'`
set    odir=`${progex} ${parafile} out_dir        | awk '{ print $2}'`



# Create the needed directories
if ( ${step} == "inst" ) then  
  if ( ! -d ${idir} ) mkdir ${idir}
  if ( ! -d ${rdir} ) mkdir ${rdir}
  if ( ! -d ${odir} ) mkdir ${odir}
endif

# Create the needed directories
if ( ${step} == "help" ) then 
  echo 
  echo "Installation"
  echo "   inst:   Creates the input, run and output directory; copy template parameter file"
  echo 
  echo "Sample case study"
  echo "   sample: Copy all files for a sample case study"
  echo 
  echo "Preparing input files [prep]"
  echo "   prep0:  Calculate S file with PV and TH [ P -> S ]"
  echo "   prep1:  Interpolate onto height levels [ P,Z -> H ]"
  echo "   prep2:  Rotate into local cartesian coordinate system [ H -> R ]"
  echo "   prep3:  Add TH,PV,NSQ and RHO to the data file [ -> R ]"
  echo "   prep4:  Define modified and anomaly PV field and boundary values [ -> R ]"
  echo "   prep5:  Reduce the domain size and split the input files [ R -> ORG,MOD,ANO,REF]" 
  echo "   prep6:  Add the reference profile [MOD -> REF]"
  echo "   prep7:  Add coastlines to REF file [-> REF]"
  echo "   prep8:  Move the files to the run directory [ORG,MOD,ANO,REF -> ]"
  echo
  echo "Perform the PV inversion [pvin]"
  echo "   pvin1:  Add NSQ, TH, RHO, and PV [ -> MOD]"
  echo "   pvin2:  Change Ertel's PV anomaly into one of quasi-geostrophic PV [MOD -> ANO]"
  echo "   pvin3:  Inversion of quasi-geostrophic PV anomaly [ANO -> ANO]"
  echo "   pvin4:  Subtract anomaly from MOD file [MOD-ANO -> MOD]"
  echo "   pvin5:  Keep iterative steps if save flag is set"
  echo
  echo "Postprocessing [post]"
  echo "   post1:  Copy needed files from input and run directory [P,MOD -> ]" 
  echo "   post2:  Rotate from quasi-cartesian coordinate frame to lat/lon system [ MOD -> GEO ]"
  echo "   post3:  Bring modified fields back to P file [ P+GEO -> P ]" 
  echo "   post4:  Calculate S file with PV and TH [ P -> S ]"
  echo "   post5:  Make clean"
  echo "   post6:  Calculate difference to original P and S file"
  echo
  echo "Diagnostic Tools"
  echo "   diag1:  Check the consitency of the boundary conditions"
  echo "   diag2:  Calculate the quasi-geostrophic PV (diag2 [ORG|MOD|ANO])"
  echo "   diag3:  Get the difference between two files (diag3 [ORG|MOD|ANO] [ORG|MOD|ANO])"
  echo "   diag4:  Calculate geopotential with hydrostatic equation (diag4 [P|ORG|MOD|ANO])"
  echo "   diag5:  Q vector analysis; geostrophic wind balance [ORG|MOD|ANO] [GEO|DIV_UV|QVEC|DIV_Q]"
  echo
  echo "Modifying Pv anomalies"
  echo "   mod1:  Stretching and amplitude changes for anomalies [MOD]" 
endif


# Copy sample files (if specified)
if ( ${step} == "sample" ) then 
  \cp ${sampledir}/P20060116_18 ${idir}
  \cp ${sampledir}/eraint_cst   ${idir}
  \cp ${sampledir}/Z20060116_18 ${idir}
  \cp ${sampledir}/pressure_cst ${idir}
endif

# ---------------------------------------------------------------------------
# Preparatory steps
# ---------------------------------------------------------------------------

# Change to data directory
cd ${idir}

# Step 0: Calculate S file with PV and TH (not in prep mode)
if ( ${step} == "prep0" ) then
    \rm -f S${date}
    /home/aswin/potentialVorticity/p2s/p2s P${date} TH PV    
    if ( $? == 0 ) then
	echo OK
    else 
	echo FAIL
    endif
endif

# Step 1: Interpolate onto height levels
if ( ${step} == "prep1" | ${step} == "prep" ) then
  echo "P${date}"                        >! fort.10
  echo "Z${date}"                        >> fort.10
  echo "H${date}"                        >> fort.10
  ${bdir}/inversion.perl ${parafile} p2z >> fort.10
    if ( $? == 0 ) then
	echo PERL SCRIPT OK
    else 
	echo FAIL
    endif

  echo "U        U       P${date}"       >> fort.10
  echo "V        V       P${date}"       >> fort.10
  echo "T        T       P${date}"       >> fort.10
  echo "Q        Q       P${date}"       >> fort.10
  ${bdir}/prep/p2z
  \rm -f fort.10 
endif

# Step 2: Rotate into local cartesian coordinate system 
if ( ${step} == "prep2" | ${step} == "prep" ) then
  echo "H${date}"                                >! fort.10
  echo "R${date}"                                >> fort.10 
  ${bdir}/inversion.perl ${parafile} rotate_grid >> fort.10
  echo "5"                                       >> fort.10
  echo "ORO"                                     >> fort.10              
  echo "U.V"                                     >> fort.10
  echo "P"                                       >> fort.10
  echo "T"                                       >> fort.10
  echo "Q"                                       >> fort.10
  ${bdir}/prep/rotate_grid
  \rm -f fort.10
endif

# Step 3: Add TH,PV,NSQ and RHO to the data file
if ( ${step} == "prep3" | ${step} == "prep" ) then
  echo "TH"         >! fort.10
  echo "R${date}"   >> fort.10
  echo "R${date}"   >> fort.10
  echo "5       "   >> fort.10
  ${bdir}/prep/z2s

  echo "PV"         >! fort.10
  echo "R${date}"   >> fort.10
  echo "R${date}"   >> fort.10
  echo "5       "   >> fort.10
  ${bdir}/prep/z2s

  echo "NSQ"        >! fort.10
  echo "R${date}"   >> fort.10
  echo "R${date}"   >> fort.10
  echo "5       "   >> fort.10
  ${bdir}/prep/z2s

  echo "RHO"        >! fort.10
  echo "R${date}"   >> fort.10
  echo "R${date}"   >> fort.10
  echo "5       "   >> fort.10
  ${bdir}/prep/z2s
  \rm -f fort.10
endif


# Step 4: Set the modified PV field and boundary values
if ( ${step} == "prep4" | ${step} == "prep" ) then
  echo "R${date}"                                >! fort.10
  ${bdir}/inversion.perl ${parafile} def_anomaly >> fort.10
  ${bdir}/prep/def_anomaly
  \rm -f fort.10 
endif

# Step 5: Reduce the domain size and split the input files
if ( ${step} == "prep5" | ${step} == "prep" ) then
  echo "PV        PV        R${date}    ORG_${date}" >! fort.10
  echo "U         U         R${date}    ORG_${date}" >> fort.10
  echo "V         V         R${date}    ORG_${date}" >> fort.10
  echo "TH        TH        R${date}    ORG_${date}" >> fort.10
  echo "Q         Q         R${date}    ORG_${date}" >> fort.10
  echo "P         P         R${date}    ORG_${date}" >> fort.10
  echo "T         T         R${date}    ORG_${date}" >> fort.10
  echo "PV_FILT   PV_AIM    R${date}    MOD_${date}" >> fort.10
  echo "U         U         R${date}    MOD_${date}" >> fort.10
  echo "V         V         R${date}    MOD_${date}" >> fort.10
  echo "Q         Q         R${date}    MOD_${date}" >> fort.10
  echo "P         P         R${date}    MOD_${date}" >> fort.10
  echo "T         T         R${date}    MOD_${date}" >> fort.10
  echo "NSQ       NSQ       R${date}    MOD_${date}" >> fort.10
  echo "RHO       RHO       R${date}    MOD_${date}" >> fort.10
  echo "TH        TH        R${date}    MOD_${date}" >> fort.10
  echo "PV_ANOM   PV        R${date}    ANO_${date}" >> fort.10
  echo "TH_ANOM   TH        R${date}    ANO_${date}" >> fort.10
  echo "UU_ANOM   U         R${date}    ANO_${date}" >> fort.10
  echo "VV_ANOM   V         R${date}    ANO_${date}" >> fort.10
  echo "ORO       ORO       R${date}    REF_${date}" >> fort.10
  echo "X         X         R${date}    REF_${date}" >> fort.10
  echo "Y         Y         R${date}    REF_${date}" >> fort.10
  echo "LAT       LAT       R${date}    REF_${date}" >> fort.10
  echo "LON       LON       R${date}    REF_${date}" >> fort.10
  echo "CORIOL    CORIOL    R${date}    REF_${date}" >> fort.10
  ${bdir}/prep/cutnetcdf
  \rm -f fort.10 
endif

# Step 6: Add the reference profile 
if ( ${step} == "prep6" | ${step} == "prep" ) then
  \rm -f fort.10 
  echo "MOD_${date}"                       >! fort.10
  echo "REF_${date}"                       >> fort.10
  ${bdir}/inversion.perl ${parafile} ref   >> fort.10
  ${bdir}/prep/ref_profile

  \rm -f fort.10 
endif

# Step 7: Add coastlines to REF file
if ( ${step} == "prep7" | ${step} == "prep" ) then
  \rm -f fort.10 
  echo \"REF_${date}\"                         >! fort.10
  echo \"${coastfile}\"                        >> fort.10
  ${bdir}/inversion.perl ${parafile} coastline >> fort.10
  ${bdir}/prep/coastline
  \rm -f fort.10 
endif

# Step 8: Move the files to the run directory
if ( ${step} == "prep8" | ${step} == "prep" ) then
  cd /home/aswin/potentialVorticity/data/inp
  \mv MOD_${date} MOD_${date}_cst ${rdir}
  \mv ORG_${date} ORG_${date}_cst ${rdir}
  \mv ANO_${date} ANO_${date}_cst ${rdir}
  \mv REF_${date} REF_${date}_cst ${rdir}
  \rm -f fort.10
endif

# ---------------------------------------------------------------------------
# Inversion
# ---------------------------------------------------------------------------

# Change to data directory
cd ${rdir}

# Start loop
set count=0
loop:

# Step 1: Add NSQ, TH, RHO, and PV to MOD file, take grid from REF file
if ( ${step} == "pvin1" | ${step} == "pvin" ) then
  echo "NSQ"           >! fort.10
  echo "MOD_${date}"   >> fort.10
  echo "REF_${date}"   >> fort.10
  echo "5          "   >> fort.10
  ${bdir}/prep/z2s
  echo "RHO"           >! fort.10
  echo "MOD_${date}"   >> fort.10
  echo "REF_${date}"   >> fort.10
  echo "5          "   >> fort.10
  ${bdir}/prep/z2s
  echo "TH"            >! fort.10
  echo "MOD_${date}"   >> fort.10
  echo "REF_${date}"   >> fort.10
  echo "5          "   >> fort.10
  ${bdir}/prep/z2s
  echo "PV"            >! fort.10
  echo "MOD_${date}"   >> fort.10
  echo "REF_${date}"   >> fort.10
  echo "5          "   >> fort.10
  ${bdir}/prep/z2s
  \rm -f fort.10
endif

# Step 2: Change Ertel's PV anomaly into an anomaly of quasi-geostrophic PV
if ( ${step} == "pvin2" | ${step} == "pvin" ) then
  echo "MOD_${date}"   >! fort.10
  echo "REF_${date}"   >> fort.10
  echo "ANO_${date}"   >> fort.10
  ${bdir}/pvin/pv_to_qgpv
  \rm -f fort.10 
endif



# Step 3: Inversion of quasi-geostrophic PV anomaly with Neumann boundary 
if ( ${step} == "pvin3" | ${step} == "pvin" ) then
  echo "ANO_${date}"   >! fort.10
  echo "REF_${date}"   >> fort.10
  ${bdir}/pvin/inv_cart
  \rm -f fort.10 
endif

# Step 4: Prepare the output of the inversion for next iteration step
if ( ${step} == "pvin4" | ${step} == "pvin" ) then
  \rm -f fort.10 
  echo "MOD_${date}"                                >! fort.10
  echo "ANO_${date}"                                >> fort.10
  ${bdir}/inversion.perl ${parafile} prep_iteration >> fort.10  
  ${bdir}/pvin/prep_iteration
  \rm -f fort.10 
endif

# Step 5: Keep iterative steps if save flag is set
if ( ${step} == "pvin5" | ${step} == "pvin" ) then
  if ( "${save}" == "yes" ) then
    set pre=''
    if ( ${count} < 10 ) then
      set pre='0'
    endif
    \cp MOD_${date} MOD_${date}_${pre}${count}
    \cp ANO_${date} ANO_${date}_${pre}${count}
  endif
endif

# End loop for iterations
if ( ${step} == "pvin" ) then
  @ count = ${count} + 1
  if ( ${count} < ${nofiter} ) goto loop
endif

# ---------------------------------------------------------------------------
# Postprocessing
# ---------------------------------------------------------------------------

# Change to output directory
cd ${odir}

# Step 1: Copy needed files from input and run directory
if ( ${step} == "post1" | ${step} == "post" ) then
    ln -sf ${rdir}/MOD_${date}       .
    ln -sf ${rdir}/MOD_${date}_cst   .
    ln -sf ${rdir}/ORG_${date}       .
    ln -sf ${rdir}/ORG_${date}_cst   .
    \cp ${idir}/P${date}             .
    set cfn = `/home/aswin/potentialVorticity/getcfn/getcfn.sh P${date}`
    \cp ${idir}/${cfn}               .
endif

# Step 2:  Rotate from quasi-cartesian coordinate frame to lat/lon system
if ( ${step} == "post2" | ${step} == "post" ) then
    echo "MOD_${date}"                                >! fort.10
    echo "GEO.MOD_${date}"                            >> fort.10   
    ${bdir}/inversion.perl ${parafile} rotate_lalo    >> fort.10
    echo "3"                                          >> fort.10
    echo "T"                                          >> fort.10     
    echo "U.V"                                        >> fort.10
    echo "P"                                          >> fort.10
    ${bdir}/post/rotate_lalo
    \rm -f fort.10 
    echo "ORG_${date}"                                >! fort.10
    echo "GEO.ORG_${date}"                            >> fort.10   
    ${bdir}/inversion.perl ${parafile} rotate_lalo    >> fort.10
    echo "3"                                          >> fort.10
    echo "T"                                          >> fort.10     
    echo "U.V"                                        >> fort.10
    echo "P"                                          >> fort.10
    ${bdir}/post/rotate_lalo
    \rm -f fort.10 
endif

# Step 3: Bring modified fields back to P file
if ( ${step} == "post3" | ${step} == "post" ) then
    \rm -f fort.10 
    echo "P${date}"                                   >! fort.10
    echo "GEO.MOD_${date}"                            >> fort.10 
    echo "GEO.ORG_${date}"                            >> fort.10  
    ${bdir}/inversion.perl ${parafile} add2p          >> fort.10
    ${bdir}/post/add2p
    \rm -f fort.10 
endif

# Step 4: Calculate S file with PV and TH
if ( ${step} == "post4" | ${step} == "post" ) then
    \rm -f S${date}
    ${bdir}/p2s/p2s P${date} TH PV    
endif

# Step 5: Make clean
if ( ${step} == "post5" | ${step} == "post" ) then
    \rm -f MOD_${date} MOD_${date}_cst
    \rm -f GEO_${date} GEO_${date}_cst
endif

# Step 6: Get difference of original and modified P,S files
if ( ${step} == "post6" | ${step} == "post" ) then
   \rm -f levs
   echo "950." >! levs
   echo "925." >> levs
   echo "900." >> levs
   echo "875." >> levs
   echo "850." >> levs
   echo "825." >> levs
   echo "800." >> levs
   echo "775." >> levs
   echo "750." >> levs
   echo "725." >> levs
   echo "700." >> levs
   echo "675." >> levs
   echo "650." >> levs
   echo "625." >> levs
   echo "600." >> levs
   echo "575." >> levs
   echo "550." >> levs
   echo "525." >> levs
   echo "500." >> levs
   echo "475." >> levs
   echo "450." >> levs
   echo "425." >> levs
   echo "400." >> levs
   echo "375." >> levs
   echo "350." >> levs
   echo "325." >> levs
   echo "300." >> levs
   echo "275." >> levs
   echo "250." >> levs
   echo "225." >> levs
   echo "200." >> levs
   echo "175," >> levs
   echo "150." >> levs
   echo "125." >> levs
   echo "100." >> levs
   echo " 75." >> levs
   echo " 50." >> levs

   \rm -f vars
   echo "U P"  >! vars
   echo "V P"  >> vars
   echo "T P"  >> vars
   echo "Z P"  >> vars
   echo "TH S" >> vars
   echo "PV S" >> vars

#   nput2p ${date} pr_cst levs vars
   ncks -A -v PS P${date} L${date}
   \mv -f L${date} L${date}.1

   \cp levs vars ${idir}/
   cd ${idir}
#   nput2p ${date} pr_cst levs vars
   ncks -A -v PS P${date} L${date}
   \rm levs
   \rm vars
   \mv L${date} ${odir}/L${date}.0
   cd ${odir} 

   ncdiff L${date}.1 L${date}.0 L${date}
   \rm -f L${date}.1
   \rm -f L${date}.0

endif


# ---------------------------------------------------------------------------
# Diagnotic Tools
# ---------------------------------------------------------------------------

# Step 1: Check the consistency of the boundary conditions (diag1)
if ( ${step} == "diag1"  ) then
    cd ${rdir}
    \rm -f fort.10 
    echo "ANO_${date}"  >! fort.10
    echo "REF_${date}"  >> fort.10   
    ${bdir}/diag/check_boundcon
    \rm -f fort.10 
endif

# Step 2: Calculate the quasi-geostrophic PV (diag2 [ORG|MOD|ANO]
if ( ${step} == "diag2"  ) then
    cd ${rdir}
    \rm -f fort.10 
    echo "${file1}_${date}"  >! fort.10
    echo "REF_${date}"       >> fort.10  
    ${bdir}/diag/calc_qgpv
    \rm -f fort.10 
endif

# Step 3: Get difference between two files (diag3 [ORG|MOD|ANO] - [ORG|MOD|ANO])
if ( ${step} == "diag3"  ) then
    cd ${rdir}
    \rm -f fort.10 
    echo "${file1}_${date}"  >! fort.10
    echo "${file2}_${date}"  >> fort.10  
    echo "DIA_${date}"       >> fort.10  
    ${bdir}/diag/difference
    \rm -f fort.10 
endif

# Step 4: Calculate geopotential with hydrostatic equation (diag4)
if ( ${step} == "diag4"  ) then
    cd ${odir}
    \rm -f fort.10 
    echo "P${date}"  >! fort.10
    echo "P${date}"  >> fort.10  
    ${bdir}/diag/hydrostatic
    \rm -f fort.10 
endif

# Step 5:  Q vector analysis / geostrophic balance check
if ( ${step} == "diag5"  ) then
    cd ${rdir}
    \rm -f fort.10
    echo "${param3}"         >! fort.10
    echo "${file1}"          >> fort.10
    echo "${file2}"          >> fort.10
    echo 1                   >> fort.10
    ${bdir}/diag/qvec_analysis
    \rm -f fort.10
endif


# ---------------------------------------------------------------------------
# Modifying tools
# ---------------------------------------------------------------------------

if ( ${step} == "mod1"  ) then
    cd ${idir}
    \rm -f fort.10 
    set commandstr=$2
    echo "R${date}"          >! fort.10  
    echo \"${commandstr}\"   >> fort.10    
    ${bdir}/spec/modify_anomaly
    \rm -f fort.10 
endif








