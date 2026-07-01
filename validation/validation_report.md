# phytreon validation (pure Python, no external tools)

## 1. Likelihood engine vs independent naive implementation (JC69)
- engine  : -103.550268
- naive   : -103.550268
- **|diff| : 8.53e-14**  (PASS, tolerance 1e-6)

## 2. NJ on additive (patristic) distances recovers the tree
- Robinson-Foulds(true, NJ) : **0**  (PASS, expected 0)

## 3. ML recovers the known A|B split
- A|B clade present : **True**  (PASS)
- logL -103.55, AIC 225.10

