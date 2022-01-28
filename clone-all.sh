#! /usr/bin/sh

## Usage: clone-all.sh will clone all of your chimes-related bitbucket repos (except for chimes_lsq)
##        into the imports directory. Prior to using this script you need to fork all
##        of the repos for chimes in your personal bitbucket space.


## BITBUCKET_USER should be your LLNL OUN.  Set in your .profile if LLNL_USER is not defined.
## BITBUCKET_SITE will vary depending on location.
CHIMES_REPOS="al_driver chimes_calculator chimes_developer_notes chimes_git_guide dlars owlqn"

## Some repos require specific branches.  Set the correct branches below.
calculator_branch="develop"

echo "Are you on an Livermore Computing system ? (y/n)"
read IS_LC

echo "Enter your LLNL OUN"
read BITBUCKET_USER

if [ "$IS_LC" == 'y' ] ; then
    BITBUCKET_SITE=ssh://git@mybitbucket.llnl.gov:7999/
else    
    BITBUCKET_SITE=https://mybitbucket.llnl.gov/scm/
fi

if [ ! -d "imports" ] ; then
	 mkdir imports
fi

cd imports

for repo in $CHIMES_REPOS ; do
    if [ ! -d $repo ] ; then
	echo "Cloning  ${BITBUCKET_SITE}~${BITBUCKET_USER}/${repo}.git"
	git clone ${BITBUCKET_SITE}~${BITBUCKET_USER}/${repo}.git $repo
	if [ $repo = "chimes_calculator" ] ; then
	    cd $repo
	    git checkout $calculator_branch
	    cd -
	fi
	if [ $? -ne 0 ] ; then
	    echo 'git clone failed'
	    exit
	fi
    else
	echo "$repo already found in imports. Not cloning"
    fi
done
cd -
