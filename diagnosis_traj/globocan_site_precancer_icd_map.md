# GLOBOCAN cancer-site precancer ICD map for UKBB trajectory modeling

Version: 2026-05-28

Purpose: define cancer-site-specific diagnosis tokens for a UKBB inpatient +
cancer-registry trajectory model: `pre-cancer -> cancer -> metastasis`, with
cancer sites aligned to Bray et al. CA: A Cancer Journal for Clinicians 2024
GLOBOCAN 2022.

This is a modeling map, not a clinical coding guideline. It is intentionally
conservative: use high-specificity precursor lesions first, then optionally add
lower-specificity risk/field-disease codes in sensitivity analyses.

## Global rules

### Cancer site classification

Use the 36 GLOBOCAN 2022 sites and ICD-10 malignant-code definitions from Bray
et al.:

| GLOBOCAN site | Cancer ICD-10 |
|---|---|
| Lip, oral cavity | C00-C06 |
| Salivary glands | C07-C08 |
| Oropharynx | C09-C10 |
| Nasopharynx | C11 |
| Hypopharynx | C12-C13 |
| Esophagus | C15 |
| Stomach | C16 |
| Colon | C18 |
| Rectum | C19-C20 |
| Anus | C21 |
| Liver, incl. intrahepatic bile ducts | C22 |
| Gallbladder | C23 |
| Pancreas | C25 |
| Larynx | C32 |
| Lung, trachea, bronchus | C33-C34 |
| Melanoma of skin | C43 |
| Non-melanoma skin cancer | C44 |
| Mesothelioma | C45 |
| Kaposi sarcoma | C46 |
| Female breast | C50 |
| Vulva | C51 |
| Vagina | C52 |
| Cervix uteri | C53 |
| Corpus uteri | C54 |
| Ovary | C56 |
| Penis | C60 |
| Prostate | C61 |
| Testis | C62 |
| Kidney incl. renal pelvis | C64-C65 |
| Bladder | C67 |
| Brain, central nervous system | C70-C72 |
| Thyroid | C73 |
| Hodgkin lymphoma | C81 |
| Non-Hodgkin lymphoma | C82-C86, C96 |
| Multiple myeloma and immunoproliferative diseases | C88, C90 |
| Leukemia | C91-C95 |

For colorectal analyses, Bray et al. combine colon, rectum, and anus as
colorectal cancer: `C18-C21`.

### ICD token tiers

- **Tier 1, direct precursor or in situ**: in situ carcinoma, intraepithelial
  neoplasia/dysplasia, named premalignant neoplasm, MGUS/smoldering myeloma.
  Use in primary models.
- **Tier 2, high-risk benign precursor where ICD cannot specify histology**:
  adenoma/polyp/cyst codes that may contain both low-risk and high-risk lesions.
  Use in sensitivity models or where histology-linked fields are available.
- **Tier 3, field disease/risk condition**: cirrhosis, chronic viral hepatitis,
  Barrett esophagus without dysplasia, atrophic gastritis, lichen sclerosus,
  immunodeficiency/HIV, etc. These are not "precancer" in the strict sequence
  sense; use as risk-context tokens, not as cancer-predecessor labels.
- **No stable ICD precursor**: do not force a precancer state unless pathology,
  screening registry, or cancer-registry morphology data are available.

UKBB inpatient ICD-10 codes are commonly stored without dots. Store both the
dotted label and no-dot token, for example `D06.9` -> `D069`.

## Site-level ICD map

| GLOBOCAN site | Cancer ICD | Primary precursor concept | ICD-10 candidate tokens | Tier | Modeling note |
|---|---:|---|---|---:|---|
| Lip, oral cavity | C00-C06 | Oral epithelial dysplasia, leukoplakia, erythroplakia, oral submucous fibrosis, carcinoma in situ | `D00.0` (`D000`), `K13.2` (`K132`), `K13.5` (`K135`), optionally `K13.7` (`K137`) | 1-2 | Use `D000` as high-specificity; `K132/K135` capture potentially malignant oral disorders but include heterogeneous severity. |
| Salivary glands | C07-C08 | No generally accepted ICD-coded precursor | none; optionally benign salivary neoplasm `D11.-` as separate benign-neoplasm history | none/2 | Do not label `D11` as precancer by default; salivary malignant transformation is uncommon and histology-dependent. |
| Oropharynx | C09-C10 | HPV-related high-grade squamous intraepithelial lesion is biologically plausible but poorly site-coded | `D00.0` (`D000`) if explicitly oral/pharyngeal CIS; HPV status codes only as risk context | 1/3 | In inpatient ICD data, true oropharyngeal precancer is usually undercaptured. |
| Nasopharynx | C11 | No stable ICD-coded precursor; EBV-associated premalignant lesions are not routinely coded | none | none | Use EBV/risk context only if laboratory data exist. |
| Hypopharynx | C12-C13 | Pharyngeal epithelial dysplasia/CIS | `D00.0` (`D000`) when pharyngeal CIS is explicit | 1 | Low sensitivity in ICD-only data. |
| Esophagus | C15 | Barrett esophagus, esophageal dysplasia/CIS | `D00.1` (`D001`), `K22.7` (`K227`), optionally achalasia `K22.0` (`K220`) as risk context | 1-3 | `D001` is direct high-grade/in situ; `K227` is best ICD token for adenocarcinoma precursor but should be stratified as Barrett context unless dysplasia is known. |
| Stomach | C16 | Gastric dysplasia/CIS, intestinal metaplasia, chronic atrophic gastritis | `D00.2` (`D002`), `K29.4` (`K294`), `K31.A*` if ICD-10-CM is present, `K31.8` (`K318`) only after local code-text validation | 1-3 | Use `D002` as strict precursor. Atrophic gastritis/metaplasia are field-risk conditions and may be poorly resolved in UK ICD-10. |
| Colon | C18 | Adenoma, advanced adenoma, serrated lesion, carcinoma in situ | `D01.0` (`D010`), `D12.0-D12.6` (`D120-D126`), `K63.5` (`K635`) | 1-2 | `D010` is high-specificity. `D12/K635` identify polyps/benign neoplasms but histology is mixed; use as Tier 2 unless pathology fields confirm adenoma/serrated dysplasia. |
| Rectum | C19-C20 | Rectal adenoma/serrated lesion/CIS | `D01.1-D01.2` (`D011-D012`), `D12.7-D12.8` (`D127-D128`), `K62.1` (`K621`) | 1-2 | Same caveat as colon. |
| Anus | C21 | Anal intraepithelial neoplasia/HSIL/CIS | `D01.3` (`D013`), optionally HPV-related dysplasia terms if local ICD extension exists | 1 | Anal HSIL/AIN may map to CIS or dysplasia codes inconsistently. |
| Liver incl. intrahepatic bile ducts | C22 | Dysplastic nodule/CIS of liver or bile duct; cirrhosis/chronic viral hepatitis as field disease | `D01.5` (`D015`) for CIS liver/biliary; `K74.*` (`K74`), `K70.3` (`K703`), `B18.*` (`B18`) as Tier 3 | 1/3 | For HCC, cirrhosis and chronic HBV/HCV are major risk-context states, not histologic precancer. Keep separate from strict precursor. |
| Gallbladder | C23 | Biliary CIS, gallbladder dysplasia; gallbladder polyp/chronic cholecystitis as risk context | `D01.5` (`D015`), `K82.4` (`K824`), optionally `K80.*` (`K80`) as risk context | 1-3 | `K824` polyps are heterogeneous; do not equate all gallstones with precancer. |
| Pancreas | C25 | PanIN, IPMN, mucinous cystic neoplasm; ICD often nonspecific | `D13.6` (`D136`), `K86.2` (`K862`) if cystic lesion, `D01.7/D01.9` if pancreatic CIS appears locally | 2 | True PanIN is rarely coded in hospital ICD. IPMN/MCN may be captured as benign/uncertain pancreatic neoplasm; validate code descriptions locally. |
| Larynx | C32 | Laryngeal dysplasia, leukoplakia/keratosis, CIS | `D02.0` (`D020`), `J38.3` (`J383`) for other vocal-cord/larynx diseases if local text supports dysplasia/leukoplakia | 1-2 | `D020` is strict; `J383` is broad and should be used only after code-text validation. |
| Lung, trachea, bronchus | C33-C34 | Bronchial squamous dysplasia/CIS; atypical adenomatous hyperplasia is not reliably ICD-coded | `D02.1-D02.2` (`D021-D022`) | 1 | High-specificity but very low sensitivity in ICD-only data. |
| Melanoma of skin | C43 | Melanoma in situ; severely dysplastic nevus is weakly coded | `D03.*` (`D03`), optionally `D22.*` (`D22`) only as benign nevus history | 1 | Use `D03` as direct pre-invasive melanoma. Do not call all nevi precancer. |
| Non-melanoma skin cancer | C44 | Actinic keratosis, Bowen disease/SCC in situ, skin CIS | `D04.*` (`D04`), `L57.0` (`L570`) | 1-2 | `D04` is high-specificity; `L570` is common and biologically relevant for SCC but not all lesions progress. |
| Mesothelioma | C45 | No accepted ICD-coded precursor | none; asbestos exposure/pleural plaque codes as exposure context if needed | none/3 | Use exposure/risk tokens separately, not precancer. |
| Kaposi sarcoma | C46 | No conventional precancer; immunosuppression/HHV-8 context | HIV/immunodeficiency codes as risk context only | 3 | Do not define a precancer sequence from ICD diagnoses. |
| Female breast | C50 | DCIS, LCIS/lobular neoplasia, atypical ductal/lobular hyperplasia | `D05.*` (`D05`), `N60.8` (`N608`) if atypical hyperplasia text is confirmed | 1-2 | `D05` includes breast CIS and is the primary precursor token. LCIS is more a risk marker than obligate precursor; keep subtype if possible. |
| Vulva | C51 | Vulvar intraepithelial neoplasia/HSIL/dVIN, vulvar CIS | `D07.1` (`D071`), `N90.0-N90.3` (`N900-N903`), lichen sclerosus `L90.0` (`L900`) as risk context | 1-3 | VIN/HSIL is direct. Lichen sclerosus is risk/field context for HPV-independent pathway. |
| Vagina | C52 | Vaginal intraepithelial neoplasia/HSIL, vaginal CIS | `D07.2` (`D072`), `N89.0-N89.3` (`N890-N893`) | 1 | Direct precursor when high-grade VaIN/HSIL. |
| Cervix uteri | C53 | CIN2/3, HSIL, adenocarcinoma in situ | `D06.*` (`D06`), `N87.1-N87.2` (`N871-N872`), optionally `N87.9` (`N879`) | 1 | Strongest ICD-coded precursor site. Use CIN2/3/HSIL and cervical CIS as strict precancer; CIN1 is lower-grade HPV lesion and should be separate. |
| Corpus uteri | C54 | Atypical endometrial hyperplasia / endometrial intraepithelial neoplasia, endometrial CIS | `D07.0` (`D070`), `N85.0` (`N850`), `N84.0` (`N840`) as polyp context | 1-2 | `N850` is useful but may include non-atypical hyperplasia; separate if 4th/5th character or text distinguishes atypia. |
| Ovary | C56 | Borderline ovarian tumor / atypical proliferative tumor; endometriosis as risk context for clear-cell/endometrioid | `D39.1` (`D391`), optionally `N80.*` (`N80`) as endometriosis context | 1/3 | `D391` is the best ICD-coded borderline/uncertain ovarian tumor token. Do not define all benign cysts as precursor. |
| Penis | C60 | Penile intraepithelial neoplasia, penile CIS | `D07.4` (`D074`), optionally `N48.0` (`N480`) balanitis xerotica obliterans as risk context | 1-3 | PeIN/CIS is direct. Chronic inflammatory penile disease is risk context only. |
| Prostate | C61 | High-grade prostatic intraepithelial neoplasia / prostate CIS | `D07.5` (`D075`), possibly `N42.3` (`N423`) dysplasia of prostate depending on local code use | 1-2 | ICD capture is weak. PIN is often pathology-only; validate before using as sequence state. |
| Testis | C62 | Germ-cell neoplasia in situ | `D07.6` (`D076`) if coded as male genital CIS/other; local morphology/pathology preferred | 1 | Rarely captured in ICD-only inpatient data. |
| Kidney incl. renal pelvis | C64-C65 | No broadly accepted ICD-coded precursor for renal-cell carcinoma; urothelial CIS for renal pelvis may appear as urinary CIS | `D09.1` (`D091`) for other urinary CIS if site text supports renal pelvis/ureter | none/1 | For RCC, avoid defining precancer from benign renal cysts. For renal pelvis urothelial carcinoma, CIS may be relevant but needs site resolution. |
| Bladder | C67 | Urothelial carcinoma in situ / flat high-grade intraurothelial lesion | `D09.0` (`D090`) | 1 | Direct and clinically meaningful precursor/non-invasive state. |
| Brain, CNS | C70-C72 | No conventional ICD-coded precursor | none | none | Do not force precancer state. |
| Thyroid | C73 | C-cell hyperplasia for medullary thyroid cancer; follicular adenoma/NIFTP not reliably mapped as precursor | `D34` (`D34`) as benign thyroid neoplasm history; `E31.2` (`E312`) MEN2 risk context | 2-3 | Thyroid precursor modeling from ICD is weak. Use RET/MEN2 genetics if available; avoid calling benign nodules universal precancer. |
| Hodgkin lymphoma | C81 | No stable ICD-coded precursor | none | none | Do not force precancer state. |
| Non-Hodgkin lymphoma | C82-C86, C96 | Monoclonal B-cell lymphocytosis for CLL/SLL; immunodeficiency/autoimmune risk context | `D72.8*` if local ICD-10-CM MBL code exists; otherwise none | none/3 | UK ICD-10 may not capture MBL. Keep risk-context tokens separate. |
| Multiple myeloma / immunoproliferative diseases | C88, C90 | MGUS, smoldering myeloma | `D47.2` (`D472`), `D47.7` (`D477`) if other lymphatic/hematopoietic uncertain behavior is used locally | 1 | MGUS is the key precursor token. Smoldering myeloma may already be within plasma-cell neoplasm coding depending on local practice. |
| Leukemia | C91-C95 | MDS/MPN/clonal hematopoiesis; ICD captures MDS/MPN better than CHIP | `D46.*` (`D46`), `D47.1` (`D471`), `D47.3` (`D473`), `D47.4` (`D474`) | 1-2 | MDS and MPN are neoplastic precursors/risk states for AML transformation, not general precursors for all leukemia. Stratify leukemia subtype if possible. |

## Recommended trajectory labels

For each participant and cancer site:

1. `site_precancer_strict`: any Tier 1 precursor code before first site-specific
   invasive cancer date.
2. `site_precancer_broad`: Tier 1 + Tier 2 before invasive cancer.
3. `site_risk_context`: Tier 3 before invasive cancer, kept as a separate token
   family.
4. `site_index_cancer`: first cancer-registry or inpatient `Cxx` event for that
   GLOBOCAN site.
5. `site_metastasis`: `C77-C79` after or on the index-cancer date.
6. `site_post_cancer`: non-metastatic diagnoses after index cancer.

Do not merge Tier 3 into the primary precancer label. For example, `K74` before
`C22` should mean `liver_risk_context_cirrhosis`, not `liver_precancer`.

## High-priority first implementation

Start with sites where ICD-coded precancer is reliable enough for trajectory
modeling:

1. Cervix: `D06`, `N871-N872`.
2. Colorectal: `D010-D013`, `D12`, `K635`, `K621`.
3. Breast: `D05`.
4. Bladder: `D090`.
5. Skin melanoma/NMSC: `D03`, `D04`, `L570`.
6. Esophagus: `D001`, `K227`.
7. Stomach: `D002`, `K294`, locally validated metaplasia/dysplasia codes.
8. Corpus uteri: `D070`, atypical `N850` if distinguishable.
9. Vulva/vagina/penis/anus: `D071-D074`, `N900-N903`, `N890-N893`.
10. Myeloma/leukemia precursors: `D472`, `D46`, selected `D47`.

## References

Cancer-site classification:

- Bray F, Laversanne M, Sung H, et al. Global cancer statistics 2022:
  GLOBOCAN estimates of incidence and mortality worldwide for 36 cancers in
  185 countries. CA: A Cancer Journal for Clinicians. 2024;74:229-263.
- WHO/IARC Global Cancer Observatory and ICD-O materials define cancer-registry
  topography using ICD-10 malignant-neoplasm site categories.

Official and high-authority precursor references used for this map:

- WHO guideline for screening and treatment of cervical pre-cancer lesions for
  cervical cancer prevention.
- NCI PDQ cervical cancer treatment: cervical precursor lesions include CIN and
  adenocarcinoma in situ.
- NCI HPV and cancer: anal, penile, vaginal, and vulvar dysplasia/intraepithelial
  neoplasia are moderate/high-grade precancerous lesions.
- NCI vulvar cancer PDQ: VIN may be a precursor to invasive vulvar squamous-cell
  carcinoma.
- NCI vaginal cancer PDQ: VaIN is classified as noninvasive squamous-cell atypia.
- NCI penile cancer PDQ: stage 0is is carcinoma in situ or PeIN.
- NCI colorectal cancer screening fact sheet: most colorectal cancers begin as
  abnormal growths/polyps; adenomas are more likely to become cancer.
- NIDDK colon polyps: removing polyps can help prevent colorectal cancer.
- NCI Barrett esophagus dictionary and esophageal cancer PDQ: Barrett esophagus
  is a premalignant condition for esophageal adenocarcinoma.
- Canadian Cancer Society stomach and pancreas precancer-condition summaries
  were used where NCI/WHO public pages were less explicit for ICD-level named
  precursors.
- NCI oral potentially malignant disorder dictionary and NCI DCEG oral cancer
  risk screening page: leukoplakia, erythroplakia, and oral submucous fibrosis
  are oral potentially malignant disorders/precursor lesions.
- NCI laryngeal cancer PDQ and NCBI MedGen/NCI concept records for laryngeal
  leukoplakia: laryngeal CIS and leukoplakia/keratosis may progress or coexist
  with invasive SCC.
- NCI bladder cancer stage and treatment PDQ: bladder stage 0is is carcinoma in
  situ, a flat tumor in the bladder lining.
- NCI breast cancer types and ACS pathology guidance: DCIS is noninvasive and
  can become invasive; LCIS increases future invasive breast cancer risk.
- NCI plasma-cell neoplasm PDQ: MGUS can become multiple myeloma or related
  hematologic malignancies.
- American Academy of Dermatology actinic keratosis overview: actinic keratosis
  and actinic cheilitis can become cutaneous squamous-cell carcinoma.

## Source URLs checked

- GLOBOCAN/CA paper local PDF: `文献资料补充/CA A Cancer J Clinicians - 2024 - Bray - Global cancer statistics 2022  GLOBOCAN estimates of incidence and mortality.pdf`
- WHO cervical pre-cancer guideline: https://www.who.int/publications/i/item/9789240030824
- NCI HPV and cancer: https://www.cancer.gov/about-cancer/causes-prevention/risk/infectious-agents/hpv-and-cancer
- NCI cervical PDQ: https://www.cancer.gov/types/cervical/hp/cervical-treatment-pdq
- NCI colorectal screening: https://www.cancer.gov/types/colorectal/screening-fact-sheet
- NIDDK colon polyps: https://www.niddk.nih.gov/health-information/digestive-diseases/colon-polyps
- NCI Barrett esophagus dictionary: https://www.cancer.gov/publications/dictionaries/cancer-terms/def/barrett-esophagus
- NCI esophageal PDQ: https://www.cancer.gov/types/esophageal/hp/esophageal-treatment-pdq
- NCI oral potentially malignant disorder dictionary: https://www.cancer.gov/publications/dictionaries/cancer-terms/def/oral-potentially-malignant-disorder
- NCI DCEG oral cancer screening: https://dceg.cancer.gov/research/cancer-types/oral-larynx-pharynx/oral-cancer-risk
- NCI vulvar PDQ: https://www.cancer.gov/types/vulvar/hp/vulvar-treatment-pdq
- NCI vaginal PDQ: https://www.cancer.gov/types/vaginal/hp/vaginal-treatment-pdq
- NCI penile PDQ: https://www.cancer.gov/types/penile/patient/penile-treatment-pdq
- NCI bladder stages: https://www.cancer.gov/types/bladder/stages
- NCI breast cancer types: https://www.cancer.gov/types/breast/breast-cancer-types
- ACS DCIS pathology guide: https://www.cancer.org/cancer/diagnosis-staging/tests/pathology-reports/breast-pathology/ductal-carcinoma-in-situ.html
- NCI plasma-cell neoplasms PDQ: https://www.cancer.gov/types/myeloma/patient/myeloma-treatment-pdq
- AAD actinic keratosis: https://www.aad.org/public/diseases/skin-cancer/actinic-keratosis-overview
- WHO ICD-O: https://www.who.int/standards/classifications/other-classifications/international-classification-of-diseases-for-oncology
- CDC ICD-10-CM official page: https://www.cdc.gov/nchs/icd/icd-10-cm/index.html
