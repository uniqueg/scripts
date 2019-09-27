#!/bin/sh

# Alexander Kanitz, Biozentrum, University of Basel
# alexander.kanitz@alumni.ethz.ch
# 19-AUG-2015

# TODO Include CLI arguments
# TODO Make quiet an option
# TODO Make repository an option
# TODO Include remote push and make an option

projectPath=${1:-"$PWD/newProject"}
repository=${2:-"ssh://git@git.scicore.unibas.ch:2222/CodeBits/directoryStructureNewProjects.git"}

gitPath=`which "git" 2> /dev/null`

if [ "$gitPath" = "" ]; then
    echo -e "[ERROR] git is required but appears not to be installed. Obtain git from 'https://git-scm.com/downloads' and make sure it is available in your \$PATH.\nExecution aborted!"
    exit 1
else
    git clone --quiet "$repository" "$projectPath"
    cd "$projectPath"
    rm -rf ".git"
    git init --quiet --template .gitTemplate
    git add .
    git commit --quiet --message "init"
fi
