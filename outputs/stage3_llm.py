import json
import requests

# Lê as proteínas de teste
with open("test_proteins.json", "r") as f:
    test_proteins = json.load(f)

# URL do Ollama
ollama_url = "http://localhost:11434/api/generate"

# Ficheiro de output onde vamos guardar os resultados
output_file = "stage3_results.json"

# Prompt
prompt = """You are not predicting protein function. You are summarizing and harmonizing provided evidence only. Do not introduce biological claims not present in the input. Do not assign EC numbers unless present in the input evidence. Do not treat STRING enrichment as direct protein-level evidence. Mention uncertainty and conflicts. Return valid JSON only.

You are a biological annotation harmonization assistant.

Task: Given the evidence packet for one protein, produce a cautious protein-level functional summary.

Rules: 1. Use only evidence present in the input. 2. Do not invent functions, pathways, EC numbers, locations, domains, or organisms.
4. Prioritize primary evidence over secondary evidence. 5. Mention weak or conflicting evidence separately. 6. If evidence is insufficient, say so. 7. Return valid JSON matching the required schema.

Input: {protein_evidence_packet}

Output schema: {{
  "protein_id": "",
  "recommended_annotation": "",
  "confidence": "strong|moderate|weak|insufficient|conflicting",
  "summary": "",
  "main_evidence": [],
  "weak_evidence": [],
  "conflicts": [],
  "warnings": [],
  "manual_review_recommended": true
}}
"""

# Lista para acumular todos os resultados
all_results = []

# Itera sobre cada proteína de teste
total = len(test_proteins)
for i, test_protein in enumerate(test_proteins, start=1):
    accession = test_protein["accession"]
    print(f"[{i}/{total}] A processar {accession}...")

    # Constrói o evidence packet
    evidence_packet = {
        "protein_id": accession,
        "protein_name": test_protein["identity"].get("protein_name"),
        "gene_name": test_protein["identity"].get("gene_name"),
        "organism": test_protein["identity"].get("organism", {}).get("scientific_name"),
        "go_annotations": test_protein.get("go_annotations", {}),
        "ec_numbers": test_protein.get("enzymatic", {}).get("ec_numbers", []),
    }

    # Substitui {protein_evidence_packet} no prompt pelo evidence packet real
    protein_prompt = prompt.format(protein_evidence_packet=json.dumps(evidence_packet))

    # Envia o evidence packet para o Ollama
    try:
        response = requests.post(ollama_url, json={
            "model": "gpt-oss",
            "prompt": protein_prompt,
            "stream": False
        }, timeout=600)

        if response.status_code == 200:
            result = response.json()
            llm_output = result.get("response", "")
            print(f"   ✓ {accession} processada com sucesso")
            all_results.append({
                "accession": accession,
                "test_group": test_protein.get("_test_group"),
                "evidence_packet": evidence_packet,
                "llm_response": llm_output,
                "status": "success"
            })
        else:
            print(f"   ✗ {accession} falhou com status {response.status_code}")
            all_results.append({
                "accession": accession,
                "test_group": test_protein.get("_test_group"),
                "status": "failed",
                "error": f"HTTP {response.status_code}"
            })

    except Exception as e:
        print(f"   ✗ {accession} falhou com erro: {e}")
        all_results.append({
            "accession": accession,
            "test_group": test_protein.get("_test_group"),
            "status": "error",
            "error": str(e)
        })

    # Guarda o ficheiro a cada iteração (para não perder progresso)
    with open(output_file, "w", encoding="utf-8") as f:
        json.dump(all_results, f, indent=2, ensure_ascii=False)

print(f"\nConcluído! Resultados guardados em {output_file}")
print(f"Total de proteínas processadas: {len(all_results)}/{total}")