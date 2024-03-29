#!/bin/bash
#=======================================================================
# mymake.sh example --Zheng Gong, 2020-03-04(yy/mm/dd)
#=======================================================================

# The below line is needed to be modified if necessary.
SRC="./src"
CompilingLog="CompilationLog.txt"

#-----------------------------------------------------------------------
# Normally no need to change anything below.
PathCurrent="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
CompilingLog=$PathCurrent/$CompilingLog
TimeString=$(date  "+%Y-%m-%d %H:%M:%S")
rm -rf $CompilingLog; touch $CompilingLog
echo                                                                  | tee -a $CompilingLog
echo "!========================*- CP3d -*========================!"   | tee -a $CompilingLog
echo "!                                                          !"   | tee -a $CompilingLog
echo "!          CP3d:    Channel-Particle 3d                    !"   | tee -a $CompilingLog
echo "!          Version: 1.0                                    !"   | tee -a $CompilingLog
echo "!          Author:  Zheng Gong                             !"   | tee -a $CompilingLog
echo "!          E-mail:  gongzheng_justin@outlook.com           !"   | tee -a $CompilingLog
echo "!                                                          !"   | tee -a $CompilingLog
echo "!====================*- Fortran 95/03 -*===================!"   | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog
echo "  Source  Path:   "$SRC                                         | tee -a $CompilingLog
echo "  Current Path:   "$PathCurrent                                 | tee -a $CompilingLog
echo "  Compiling Time: "$TimeString                                  | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog

# Set exe name
echo "  Which EXE do you want to compile? "                           | tee -a $CompilingLog
echo "     1: channel2nd"                                             | tee -a $CompilingLog
echo "     2: channel4th"                                             | tee -a $CompilingLog
if [[ -n $1 ]]; then
  strTemp=$1
  EXE=${strTemp:5}
else
  read -p "  Please type a EXE index(1 or 2): " id_exe
  echo    "  Please type a EXE index(1 or 2): "$id_exe >> $CompilingLog
  if [ "$id_exe" == 1 ]; then
    EXE="channel2nd"
  elif [ "$id_exe" == 2 ]; then
    EXE="channel4th"
  else
    echo "  Sorry, EXE type cannot be recognized."                    | tee -a $CompilingLog
    echo "  Compiling filed"                                          | tee -a $CompilingLog
    exit 1
  fi
fi
echo "  "$EXE"  will be compiled"                                     | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog

# Set compiler
echo "  Which compiler do you use? "                                  | tee -a $CompilingLog
echo "     1: Intel MPI (mpiifort)"                                   | tee -a $CompilingLog
echo "     2: gcc MPI   (mpif90). Default"                            | tee -a $CompilingLog
if [[ -n $2 ]]; then
  strTemp=$2
  CMP=${strTemp:5}
else
  read -p "  Please type a compiler index (1 or 2): " id_cmp
  echo    "  Please type a compiler index (1 or 2): "$id_cmp >> $CompilingLog
  if [ "$id_cmp" == 1 ]; then
    CMP="intel_MPI"
  else
    CMP="gcc_MPI"
  fi
fi
echo "  "$CMP"  will be used"                                         | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog

# Compile ThirdParty or not
echo "  Do you want to recompile ThirdParty? "                        | tee -a $CompilingLog
echo "     0: No need, use the old ThirdParty compilation. Default"   | tee -a $CompilingLog
echo "     1: Yes, recompile ThirdParty (Recommended for first use)"  | tee -a $CompilingLog
if [[ -n $3 ]]; then
  strTemp=$3
  Third_flag=${strTemp:0-1}
  echo    "  Please type a choice (0 or 1): "$Third_flag              | tee -a $CompilingLog
else
  read -p "  Please type a choice (0 or 1): " Third_flag
  echo    "  Please type a choice (0 or 1): "$Third_flag >> $CompilingLog
fi
if [ "$Third_flag" == 1 ]; then
  echo                                                                | tee -a $CompilingLog
  cd $SRC/ThirdParty
  echo "Compiling ThirdParty begins."                                 | tee -a $CompilingLog
  chmod a+x ./install_thirdParty.sh
  ./install_thirdParty.sh
  echo "Compiling ThirdParty done !"                                  | tee -a $CompilingLog
  echo                                                                | tee -a $CompilingLog
  cd ../..
else
  echo "  Choose to use the old ThirdParty compilation"               | tee -a $CompilingLog
fi
echo                                                                  | tee -a $CompilingLog

# Delete temporary compiling files or not
echo "  Do you want to delete temporary compiling files? "            | tee -a $CompilingLog
echo "     0: No, save them. "                                        | tee -a $CompilingLog
echo "     1: Yes,delete them. Default"                               | tee -a $CompilingLog
if [[ -n $4 ]]; then
  strTemp=$4
  DeleteFlag=${strTemp:0-1}
  echo    "  Please type a choice (0 or 1): "$DeleteFlag              | tee -a $CompilingLog
else
  read -p "  Please type a choice (0 or 1): " DeleteFlag
  echo    "  Please type a choice (0 or 1): "$DeleteFlag >> $CompilingLog
fi
if [ "$DeleteFlag" != 0 ]; then
  echo "  Choose to DELETE temporary compiling files"                 | tee -a $CompilingLog
else
  echo "  Choose to SAVE temporary compiling files"                   | tee -a $CompilingLog
fi
echo                                                                  | tee -a $CompilingLog

# Compile begins
echo  "!==================*- Compiling begins -*=================!"   | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog
rm -fr $EXE
cd $SRC
if [ "$DeleteFlag" != 0 ]; then
  make -f "make_"$EXE clean >&/dev/null
fi
if [ "$EXE" == "channel2nd" ]; then
  CFD_DEFS_Add=""
  if [[ -n $5 ]]; then
    strTemp=$5
    CFD_DEFS_Add=${strTemp:13}
  fi  
  make -f "make_"$EXE CMP=$CMP exeName=$EXE CFD_DEFS_Add=$CFD_DEFS_Add 2>&1 | tee -a $CompilingLog
elif [ "$EXE" == "channel4th" ]; then
  CFD_DEFS_Add=""
  IsSolveScalar=""
  if [[ -n $5 ]]; then
    strTemp=$5
    CFD_DEFS_Add=${strTemp:13}
  fi  
  if [[ -n $6 ]]; then
    strTemp=$6
    IsSolveScalar=${strTemp:15}
  fi
  if [ "$IsSolveScalar" == "" ]; then
    make -f "make_"$EXE CMP=$CMP exeName=$EXE CFD_DEFS_Add=$CFD_DEFS_Add 2>&1 | tee -a $CompilingLog  
  else
    make -f "make_"$EXE CMP=$CMP exeName=$EXE CFD_DEFS_Add=$CFD_DEFS_Add IsSolveScalar=$IsSolveScalar 2>&1 | tee -a $CompilingLog  
  fi
else
  echo  $EXE" wrong, please check !!!"                                | tee -a $CompilingLog
fi
echo                                                                  | tee -a $CompilingLog
mv $EXE $PathCurrent                              
if [ $? -ne 0 ]; then
  if [ "$DeleteFlag" != 0 ]; then
    make -f "make_"$EXE clean >&/dev/null
  fi
  echo  $EXE" CANNOT be compiled correctly, please check !!!"         | tee -a $CompilingLog
else
  if [ "$DeleteFlag" != 0 ]; then
    make -f "make_"$EXE clean >&/dev/null
  fi
  echo  $EXE" has been compiled normally. Enjoy !!!"                  | tee -a $CompilingLog
  cd ..
  chmod a+x ./$EXE
fi
echo                                                                  | tee -a $CompilingLog
echo  "!===================*- Compiling ends -*==================!"   | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog
