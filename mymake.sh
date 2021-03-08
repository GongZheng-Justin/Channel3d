#!/bin/bash
#=======================================================================
# mymake.sh example --Zheng Gong, 2020-03-04(yy/mm/dd)
#=======================================================================

# The below line is needed to be modified if necessary.
SRC="./src/"
CompilingLog="Channel3d_CompilingLog.txt"

#-----------------------------------------------------------------------
# Normally no need to change anything below
PathCurrent="$(dirname "$(readlink -f ${BASH_SOURCE[0]})")"
CompilingLog=$PathCurrent/$CompilingLog
TimeString=$(date  "+%Y-%m-%d %H:%M:%S")
rm -rf $CompilingLog; touch $CompilingLog
echo                                                                | tee -a $CompilingLog
echo "!======================*- Channel3d -*=====================!" | tee -a $CompilingLog
echo "!                                                          !" | tee -a $CompilingLog
echo "!          CP3d:    Channel3d                              !" | tee -a $CompilingLog
echo "!          Version: 1.0                                    !" | tee -a $CompilingLog
echo "!          Author:  Zheng Gong                             !" | tee -a $CompilingLog
echo "!          E-mail:  gongzheng_justin@outlook.com           !" | tee -a $CompilingLog
echo "!                                                          !" | tee -a $CompilingLog
echo "!====================*- Fortran 95/03 -*===================!" | tee -a $CompilingLog
echo                                                                | tee -a $CompilingLog
echo "  Source  Path:   "$SRC                                       | tee -a $CompilingLog
echo "  Current Path:   "$PathCurrent                               | tee -a $CompilingLog
echo "  Compiling Time: "$TimeString                                | tee -a $CompilingLog
echo                                                                | tee -a $CompilingLog
echo "  Which compiler do you use? "                                | tee -a $CompilingLog
echo "     1: Intel MPI ( mpiifort )"                               | tee -a $CompilingLog
echo "     2: gcc MPI   ( mpif90   )"                               | tee -a $CompilingLog
read -p "  please type a compiler index (1 or 2): " id_cmp
echo    "  please type a compiler index (1 or 2): "$id_cmp >> $CompilingLog

if [ $id_cmp -eq 1 ]; then
  CMP="intel_MPI"
elif [ $id_cmp -eq 2 ]; then
  CMP="gcc_MPI"
else
  echo "  Sorry, compiler type cannot be recognized."               | tee -a $CompilingLog
  echo "  Compiling filed"                                          | tee -a $CompilingLog
  exit 1
fi
echo "  "$CMP"  will be used"                                       | tee -a $CompilingLog
echo                                                                | tee -a $CompilingLog

echo "  Which EXE do you want to compile? "                         | tee -a $CompilingLog
echo "     1: channel2nd"                                           | tee -a $CompilingLog
echo "     2: channel4th"                                           | tee -a $CompilingLog
read -p "  please type a EXE index(1 or 2): " id_exe
echo    "  please type a EXE index(1 or 2): "$id_exe >> $CompilingLog

if [ $id_exe -eq 1 ]; then
  EXE="channel2nd"
elif [ $id_exe -eq 2 ]; then
  EXE="channel4th"
else
  echo "  Sorry, EXE type cannot be recognized."                      | tee -a $CompilingLog
  echo "  Compiling filed"                                            | tee -a $CompilingLog
  exit 2
fi
echo "  "$EXE"  will be compiled"                                     | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog
echo  "!==================*- Compiling begins -*=================!"   | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog
rm -fr $EXE
cd $SRC
make -f "make_"$EXE clean >&/dev/null 
make -f "make_"$EXE CMP=$CMP   2>&1                                   | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog
mv $EXE $PathCurrent                              
if [ $? -ne 0 ]; then
  make -f "make_"$EXE clean >&/dev/null
  echo  $EXE" CANNOT be compiled correctly, please check !!!"         | tee -a $CompilingLog
else
  make -f "make_"$EXE clean >&/dev/null
  echo  $EXE" has been compiled normally. Enjoy !!!"                  | tee -a $CompilingLog
  cd ..
  chmod a+x ./$EXE
fi
echo                                                                  | tee -a $CompilingLog
echo  "!===================*- Compiling ends -*==================!"   | tee -a $CompilingLog
echo                                                                  | tee -a $CompilingLog
