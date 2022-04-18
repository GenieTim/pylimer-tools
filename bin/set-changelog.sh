#!/bin/sh
# replace a placeholder for next v content with next v content

# setup variables
SCRIPT_DIR=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
PY_VERSION=$(sed -nr 's/[ \t]*version="([.0-9\.]*)",[ \t]*/\1/p' setup.py)

# synchronize the version also to the CMake file
sed -ri 's/project\(pylimer_tools VERSION "[0-9.]*"\)/project(pylimer_tools VERSION "'"$PY_VERSION"'")/g' CMakeLists.txt

# construct the changelog since the last change
CHANGES=$(npx -y changelog-maker --commit-url="https://github.com/BernhardWebstudio/DataShot_DesktopApp/commit/{ref}")
# switch to directory with changelog file
cd "$SCRIPT_DIR/../" || exit
EXISTING_CHANGES=$(cat CHANGELOG.md)
EXISTING_CHANGES=$(echo "$EXISTING_CHANGES" | sed "s/# Changelog//g")
printf "# Changelog\n\n## v%s\n%s\n%s" "$PY_VERSION" "$CHANGES" "$EXISTING_CHANGES" >CHANGELOG.md
