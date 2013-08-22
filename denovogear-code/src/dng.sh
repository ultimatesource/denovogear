#!/usr/bin/sh

# version info
version_num=0.1

# set default libexec location
# under certail conditions this may fail to work, so an installer should set
# this to the proper installation directory, e.g. $prefix/libexec
: ${DNG_LIBEXEC_DIR:="$( cd "$( dirname "$0" )" && pwd )"/../libexec}

# Internal commands
int_cmds="help version citation"

submatch () { case "$2" in *$1*) return 0 ;; *) return 1 ;; esac ; }

function check_libexec() {
    if [ ! -d $DNG_LIBEXEC_DIR ]; then
	echo "ERROR: DNG_LIBEXEC_DIR ($DNG_LIBEXEC_DIR) does not exist."
	echo "Solution: Try setting it in your environment."
	exit 1
    fi
}

# calls external command
function call_command_ex() {
    check_libexec
    cmd=$1
    shift
    cmdpath="$DNG_LIBEXEC_DIR/dng-$cmd"
    if [ ! -x $cmdpath  ]; then
	echo "ERROR: command '$cmd' not recognized"
	return 1
    fi
    $cmdpath $@
    return $?
}

# calls commands
function call_command() {
    cmd=$1
	# check if the command is internal
	# other wise run external command
    if submatch "$cmd" "$int_cmds"; then
	dng_$cmd $@
    else
	call_command_ex $@
    fi
    return $?
}

function dng_version() {
    echo "DeNovoGear Version $version_num"
    return 0
}

function dng_help() {
    check_libexec
    
    if [ -n "$2" ]; then
	if [ "$2" = "help" ]; then
	    echo "USAGE: dng help"
	    echo "       dng help command"
	    echo ""
	    echo "Description: Displays help information."
	    echo "    The first form displays a list of available commands."
	    echo "    The second form displays the help information for a command."
	else
	    call_command $2 "help"
	fi
	return $?
    fi
    
    ext_cmds=`find $DNG_LIBEXEC_DIR -maxdepth 1 -type f -perm +111 -name "dng-*" -print | sed -e 's/.*dng-//'`
    cmds="$int_cmds $ext_cmds"
    sorted=`echo $cmds | tr ' ' '\n' | sort | tr '\n' ' '`
    
    echo "The following commands are supported by this installation of dng:"
    
    for i in $sorted; do
	echo "    $i"
    done
    echo "For more information use 'dng help command'"
    return 0
}

function dng_usage() {
    echo "USAGE: dng command [options]"
    echo "       dng help"
    echo "       dng help command"
}

function dng_citation() {
    echo "WARNING: Not yet implemented."
    return 0
}

# print usage if no arguments were given
if [ "$#" -eq 0 ]; then
# print help message
    dng_usage
    exit 0
fi

call_command $@

exit $?

