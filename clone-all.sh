#!/bin/bash

## Usage: clone-all.sh will clone all of your chimes-related bitbucket repos (except for chimes_lsq)
##        into the imports directory. Prior to using this script you need to fork all
##        of the repos for chimes in your personal bitbucket space.
##
## BITBUCKET_USER should be your LLNL OUN.  Set in your .profile if LLNL_USER is not defined.
## BITBUCKET_SITE will vary depending on location.



FROM_INSTALL=${1-0} # If called from install.sh, don't prompt the user for LC status
IS_LC=0             # If called from install.sh, is set correctly. Otherwise, user will reset below



echo "Do your Bitbucket ChIMES forks use the default name? (y/n)"
echo "Eg.: Default fork name for chimes_calculator is  chimes_calculator-fork"

read IS_DEFAULT

if [ "$IS_DEFAULT" == 'n' ] ; then
	echo "Please modify \"REPO_FORK_NAME\" values with your fork names"
	echo "e.g.: chimes_calculator --> chimes_calculator-fork"
	echo "Then rerun this script with option 'y'"
	exit 0
fi

# Only download required packages

CHIMES_REPOS[0]=chimes_calculator	; REPO_FORK_NAME[0]=${CHIMES_REPOS[0]}-fork; TARGET_BRANCH[0]=master

# Non-neccesary packages
#CHIMES_REPOS[1]=dlars			; REPO_FORK_NAME[1]=${CHIMES_REPOS[1]}-fork; TARGET_BRANCH[1]=master
#CHIMES_REPOS[2]=chimes_developer_notes	; REPO_FORK_NAME[2]=${CHIMES_REPOS[2]}-fork; TARGET_BRANCH[2]=master
#CHIMES_REPOS[3]=chimes_git_guide	; REPO_FORK_NAME[3]=${CHIMES_REPOS[3]}-fork; TARGET_BRANCH[3]=master
#CHIMES_REPOS[4]=al_driver		; REPO_FORK_NAME[4]=${CHIMES_REPOS[4]}; TARGET_BRANCH[4]=master
#CHIMES_REPOS[5]=owlqn			; REPO_FORK_NAME[5]=${CHIMES_REPOS[5]}; TARGET_BRANCH[5]=master

if [ $FROM_INSTALL -eq 1 ] ; then
	IS_LC=$2
else

	echo "Are you on an Livermore Computing system ? (y/n)"
	read IS_LC
fi

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

echo "===="

NREPOS=${#CHIMES_REPOS[*]}

for (( repo=0; repo<${NREPOS}; repo++ )); do

    echo "... Cloning and installing... ${REPO_FORK_NAME[$repo]}"

    if [ ! -d ${CHIMES_REPOS[$repo]} ] ; then
	echo "Cloning  ${BITBUCKET_SITE}~${BITBUCKET_USER}/${REPO_FORK_NAME[$repo]}.git to folder named: ${CHIMES_REPOS[$repo]}"
	git clone ${BITBUCKET_SITE}~${BITBUCKET_USER}/${REPO_FORK_NAME[$repo]}.git ${CHIMES_REPOS[$repo]}
	    cd ${CHIMES_REPOS[$repo]}
	    git checkout ${TARGET_BRANCH[$repo]}
	    ./install.sh 
	    cd - 1&>/dev/null
	if [ $? -ne 0 ] ; then
	    echo 'git clone failed'
	    exit
	fi
    else
	echo "${REPO_FORK_NAME[$repo]} already found in imports. Not cloning."
    fi
    
    echo " ===="
done
cd - 1&>/dev/null
