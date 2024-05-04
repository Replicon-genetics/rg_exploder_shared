# Shell script to switch which Biopython_fix version is used; takes 1 parameter, defaults to current
#rootapplicationdir="/Users/caryodonnell/Documents/repositories/snowlizardz/rg_exploder/"
rootapplicationdir="/Users/caryodonnell/Documents/repositories/rg_exploder_shared/"
inputmetro="metro"
inputretro="retro"
fixhead="Biopython_fix"
under="_"
fixtail=".py"
cd $rootapplicationdir
/bin/rm Biopython_fix.py # Don't worry about 'not found' warnings
echo $1
if [ $1 == $inputretro ]
then
    fix_mid=$inputretro
else
    fix_mid=$inputmetro
fi

ln -s $fixhead$under$fix_mid$fixtail $fixhead$fixtail
