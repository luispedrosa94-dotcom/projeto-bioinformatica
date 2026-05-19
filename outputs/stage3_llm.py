import json
import requests

# Lê as proteínas de teste
with open("test_proteins.json", "r") as f:
    test_proteins = json.load(f)

# URL do Ollama
ollama_url = "http://localhost:8000/v1/complete"

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

# Itera sobre cada proteína de teste
for test_protein in test_proteins:
    # Constrói o evidence packet
    evidence_packet = {
        "protein_id": test_protein["accession"],
        "protein_name": test_protein["identity"]["protein_name"],
        "gene_name": test_protein["identity"]["gene_name"],
        "organism": test_protein["identity"]["organism"]["scientific_name"],
        "go_terms": test_protein["go_terms"],
        "ec_numbers": test_protein["ec_numbers"],
        # Adicione outros campos relevantes aqui
    }

    # Substitui {protein_evidence_packet} no prompt pelo evidence packet real
    protein_prompt = prompt.format(protein_evidence_packet=json.dumps(evidence_packet))

    # Envia o evidence packet para o Ollama
    response = requests.post(ollama_url, json={
        "prompt": protein_prompt,
    })

    # Processa a resposta do Ollama
    if response.status_code == 200:
        result = response.json()
        print(f"Resultado para a proteína {test_protein['accession']}:")
        print(result)
        print("---")
    else:
        print(f"Request failed with status {response.status_code} for protein {test_protein['accession']}")