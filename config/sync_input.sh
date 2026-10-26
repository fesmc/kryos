#!/bin/sh
# The synchronised copies of input/ (docs/running.md): the library defaults and
# tables a run reads, each a copy of a source file in a dependency checkout (or
# in yelmox), listed in input/sync.txt as <copy> <root variable> <source>.
#
#     sh config/sync_input.sh check <manifest> <ROOT>=<dir> ...   (make check)
#     sh config/sync_input.sh sync  <manifest> <ROOT>=<dir> ...   (make sync-input)
#
# check: exit 1 naming every copy that differs from its source or is missing.
# sync: copy each source over its copy (only the copies that differ).
# Both exit 1 if a source is missing or a manifest root is not given.

mode=$1
manifest=$2
shift 2

case "$mode" in
    check|sync) ;;
    *) echo "usage: sh config/sync_input.sh check|sync <manifest> <ROOT>=<dir> ..." >&2; exit 2 ;;
esac

dir=$(dirname "$manifest")

# The roots, as root_<NAME>=<dir>.
for pair in "$@"; do
    name=${pair%%=*}
    case "$name" in
        *[!A-Za-z0-9_]*|"") echo "sync_input.sh: bad root argument $pair" >&2; exit 2 ;;
    esac
    eval "root_$name=\${pair#*=}"
done

stale=""
missing=""
n=0
while read -r copy root src rest; do
    case "$copy" in ""|\#*) continue ;; esac
    eval "base=\${root_$root-}"
    if [ -z "$base" ]; then
        echo "sync_input.sh: $manifest: root $root of $copy is not given" >&2
        exit 1
    fi
    from="$base/$src"
    if [ ! -f "$from" ]; then
        missing="$missing
    input/$copy <- $from ($root = $base)"
        continue
    fi
    n=$((n + 1))
    cmp -s "$from" "$dir/$copy" && continue
    if [ "$mode" = sync ]; then
        cp "$from" "$dir/$copy" || exit 1
        echo "    input/$copy <- $from"
    else
        what="differs from $from"
        [ -f "$dir/$copy" ] || what="missing (source $from)"
        stale="$stale
    input/$copy: $what"
    fi
done < "$manifest"

if [ -n "$missing" ]; then
    echo ""
    echo "*** ERROR: source files missing (set the root on the make command line, e.g. YELMOXROOT=../yelmox):$missing"
    echo ""
    exit 1
fi
if [ -n "$stale" ]; then
    echo ""
    echo "*** ERROR: stale copies in input/ (input/sync.txt):$stale"
    echo "    Update them with: make sync-input"
    echo ""
    exit 1
fi
if [ "$mode" = check ]; then
    echo "    input/: $n synchronised copies match their sources ($manifest)."
fi
