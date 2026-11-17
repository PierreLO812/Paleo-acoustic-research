# T3 cycle — Paleo Grokbot × Antigravity

Repo local:
`C:\Users\balad\OneDrive\Documents\Paléo-acoustique\paleoacoustique`
GitHub: `https://github.com/PierreLO812/Paleo-acoustic-research`

## Loop
1. **HYPOTHESIS (Grokbot / local)**  
   `python -m pipeline.orchestrator --taxon parasaurolophus_ram14000 --mode agents`  
   → lit le `prompt_path` dans le JSON.

2. **BUILD_AND_SIM (Antigravity)**  
   `delegate_task` provider=`antigravity`, runtimeMode=`auto-accept-edits`  
   Task = contenu de `output/<run_id>/antigravity_prompt.md`  
   Attendre la fin.

3. **PEER_REVIEW (Grokbot / local)**  
   `python -m pipeline.orchestrator --taxon parasaurolophus_ram14000 --run-id <RUN_ID> --resume-review --attempt 1`

4. **RETRY_LOOP**  
   Si `AWAITING_ANTIGRAVITY_RETRY`, renvoyer le nouveau prompt à Antigravity.  
   Stop à 3 tentatives ou `VALIDATED`.

## Local smoke (sans agents)
```powershell
python -m pipeline.orchestrator --taxon parasaurolophus_ram14000 --mode local
```
