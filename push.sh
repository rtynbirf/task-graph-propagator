#!/usr/bin/env bash
# ============================================================================
#  push.sh — интеграция в GitHub БЕЗ передачи токена боту/логам.
#  Токен берётся ТОЛЬКО из окружения и НЕ остаётся в .git/config.
#
#  Использование (пользователь выполняет сам, на своей машине):
#     export GH_TOKEN=<твой новый fine-grained токен>
#     ./push.sh https://github.com/USER/REPO.git
# ============================================================================
set -euo pipefail

: "${GH_TOKEN:?Сначала: export GH_TOKEN=<токен> — токен не вводится в код и не хранится}"
REMOTE="${1:?Использование: ./push.sh <repo-url>  (например, https://github.com/USER/REPO.git)}"

git remote remove origin 2>/dev/null || true
git remote add origin "https://x-access-token:${GH_TOKEN}@${REMOTE#https://}"
git push -u origin HEAD
git remote remove origin   # токен не остаётся в .git/config

echo "OK: запушено в ${REMOTE}; токен из конфига удалён."
