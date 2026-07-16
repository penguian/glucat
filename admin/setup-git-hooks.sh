#!/bin/bash
#
# setup-git-hooks.sh: Setup Git pre-push hook for GluCat verification.
#
#    copyright            : (C) 2026 by Paul C. Leopardi
#
#    This library is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published
#    by the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with this library.  If not, see <http://www.gnu.org/licenses/>.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

PRE_COMMIT_HOOK="${ROOT_DIR}/.git/hooks/pre-commit"
PRE_PUSH_HOOK="${ROOT_DIR}/.git/hooks/pre-push"

if [ ! -d "${ROOT_DIR}/.git" ]; then
    echo "Error: ${ROOT_DIR} is not a git repository."
    exit 1
fi

# Pre-commit hook: Ultra-fast 39ms checks (core build, license headers, Ruff linter)
cat << 'EOF' > "${PRE_COMMIT_HOOK}"
#!/bin/bash
echo "=== Pre-commit fast lint gate starting ==="
python3 verify_all.py -q --fast
if [ $? -ne 0 ]; then
    echo "=== Pre-commit lint gate FAILED! Commit aborted. ==="
    exit 1
fi
echo "=== Pre-commit lint gate PASSED! Proceeding with commit. ==="
EOF

# Pre-push hook: Full verification gate (~6.6s)
cat << 'EOF' > "${PRE_PUSH_HOOK}"
#!/bin/bash
echo "=== Pre-push verification gate starting ==="
python3 verify_all.py -q --python
if [ $? -ne 0 ]; then
    echo "=== Pre-push verification gate FAILED! Push aborted. ==="
    exit 1
fi
echo "=== Pre-push verification gate PASSED! Proceeding with push. ==="
EOF

chmod +x "${PRE_COMMIT_HOOK}" "${PRE_PUSH_HOOK}"
echo "Successfully installed git hooks:"
echo "  Pre-commit: ${PRE_COMMIT_HOOK}"
echo "  Pre-push:   ${PRE_PUSH_HOOK}"
