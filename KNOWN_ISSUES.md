# KNOWN_ISSUES · task-graph-propagator · Волна ×5 (2026-11-09)

## Красные CI (было → после патча)

| Run | Что было | Диагноз | После |
|---|---|---|---|
| 42312192273 (11-06 23:29) | guardian-tripwire failure, **0 jobs**, имя рана = путь к файлу | общий корень кластера ×5: пин-склейка `# vXwith:` в tripwire.yml SHA ea5570bf (PyYAML-репро, строки 39–41) | волна ×4 переписала (385c9a5): зелёный **42516753275** @ 385c9a5c; ×5-D диспатч: run ****42527376745** (workflow_dispatch → success, 21:52, sha 385c9a5c)** |
| 42312191171 (11-06 23:29) | CI failure, **0 jobs** | тот же шаблонный косяк в ci.yml того SHA: `uses: actions/setup-python@… # v5with:` (строка 21–22) | волна ×4 переписала + добавила property-тесты (385c9a5): зелёный **42516753093** @ 385c9a5c; ×5-D диспатч: run ****42527380947** (workflow_dispatch → success, 21:52, sha 385c9a5c)** |
