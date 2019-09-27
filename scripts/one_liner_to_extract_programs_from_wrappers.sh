# Alexander Kanitz|Biozentrum|University of Basel
# 17-AUG-2015

bashBuiltins=":|.|alias|bg|bind|break|builtin|caller|cd|command|compgen|complete|compopt|continue|declare|dirs|disown|echo|enable|eval|exec|exit|export|false|fc|fg|getopts|hash|help|history|jobs|kill|let|local|logout|mapfile|popd|printf|pushd|pwd|read|readonly|return|set|shift|shopt|source|suspend|test|times|trap|true|type|typeset|ulimit|umask|unalias|unset|wait"
bashKeywords=`compgen -k | tr '\n' '|'`
#bashKeywords="next|continue|then|elif|else|fi|do|done"

sed 's/[#].*$//' "$@" | grep --only-matching --perl-regexp '(^\s*[\w\d_\-.=\[\]]+)|(;\s*[\w\d_.\-=\[\]]+)|(\|\s*[\w\d_.\-=\[\]]+)' | sed 's/^\(;\||\)\s*//' | sed 's/^\s\+//' | grep --invert-match --extended-regexp '=|[|]' | grep --invert-match --extended-regexp --line-regexp "$bashKeywords|$bashBuiltins" | sort -u
# remove everything after comment | grep first word after whitespace or semicolon | remove semicolon
# + whitespace | remove leading whitespace | remove variable assignments | remove bash keywords | 
# sort and remove duplicates

# Known issues

