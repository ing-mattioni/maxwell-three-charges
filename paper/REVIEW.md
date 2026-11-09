# Revisione scientifica e riproducibilita

Revisione interna del 30-31 ottobre 2026, assistita da OpenAI Codex.
Non costituisce una peer review esterna o una formalizzazione della dimostrazione.
Base Git esaminata: `d75df2b` (hash dopo la riscrittura della cronologia dell'8 novembre 2026,
che ha cambiato solo l'email dell'autore); le modifiche descritte sono locali.

## Valutazione

**Il lavoro merita un paper scientifico.** La domanda e precisa e il risultato
rivendicato e sostanziale: al massimo quattro equilibri per tre cariche
coulombiane positive di intensita arbitrarie. La riduzione alla geometria della
caustica, seguita dalla verifica globale delle disuguaglianze, offre una proposta
di dimostrazione assistita dal calcolatore verificabile e riproducibile.

Il manoscritto originale e stato letto sia come sorgente sia come PDF, con
ispezione visiva di tutte le sue 14 pagine. Contiene gia argomenti piu completi
delle vecchie note: blow-up dei vertici per il grado di Z, controllo della
rotazione mediante un triangolo base, confronto di Sturm con trattamento degli
estremi, e passaggio ai parametri degeneri mediante persistenza delle selle.
Questi passaggi sono stati esaminati; restano da sottoporre a un esperto esterno.

La novita va attribuita con precisione. Il limite di tre tangenti segue dal
Teorema 2 di Martinez-Maure (2008), applicato alla funzione di supporto
antiperiodica: [articolo e DOI](https://doi.org/10.1215/ijm/1254403726).
Il paper ora lo cita. Il contributo proposto riguarda l'applicazione alla mappa
delle cariche, le formule e le disuguaglianze verificate sul dominio completo;
il confronto di Sturm non va presentato come un risultato nuovo in generale.

Il confronto bibliografico include il limite sei di
[Gabrielov et al. (2026)](https://arxiv.org/abs/2607.28785), il controesempio
a cinque cariche di [Arathoon et al. (2026)](https://arxiv.org/abs/2607.27197)
e il risultato di finitezza per cariche di un segno di
[Enciso e Peralta-Salas (2026)](https://arxiv.org/abs/2609.03965).
La ricerca mirata non e una garanzia di priorita o una rassegna esaustiva.

## Correzioni integrate

| Punto | Intervento e verifica |
|---|---|
| Copertura del simplesso | `domains.py` racchiude gli estremi razionali della griglia con arrotondamento verso l'esterno. Eliminate le sottili lacune diagonali della vecchia griglia n0=12. Test con frazioni esatte. |
| Resti di Taylor | Le potenze nei resti sono ora prodotti a intervalli; eliminata la dipendenza dal calcolo non verificato delle potenze floating point. |
| Riduzione angolare | Il prodotto con multipli negativi di pi mantiene l'ordine degli estremi. Aggiunto un caso di regressione specifico. |
| Risultati incompleti | Entrambi i certificatori restituiscono exit code 1 se restano box irrisolte. Due esecuzioni con zero iterazioni verificano il comportamento. |
| Tracciabilita | Ogni JSON registra parametri, piattaforma, versioni e hash SHA256. `reproduce.py` conserva i log, controlla i conteggi e verifica che i sorgenti non cambino durante l'esecuzione. |
| Esposizione matematica | Distinti i minimi non degeneri dai casi degeneri; corretto il segno della componente trasversa nel caso collineare e il segno della traccia nella discussione delle cariche miste. |
| Portata delle verifiche | Chiarito il ruolo locale di W e del triangolo base; rimosse affermazioni eccessive sulla regolarita agli angoli singolari e sui conteggi numerici di Newton. |
| Bibliografia | Inseriti Martinez-Maure 2006/2008 e il risultato di finitezza del 2026. |
| Documento | Aggiornati risultati e tempi, spezzate formule troppo lunghe, aggiunte risorse KaTeX locali e controlli di build. |

## Esecuzioni sulla versione corretta

CPython 3.11.4 x64, NumPy 2.4.6, mpmath 1.3.0, Windows; 16 processi per
le certificazioni principali. Il runtime recuperato localmente permette questi
test, ma e incompleto: per una riproduzione esterna usare un'installazione Python
completa. Nessun pacchetto simulato e stato utilizzato.

| Esecuzione | Box elaborate | Irrisolte | Tempo del certificatore |
|---|---:|---:|---:|
| Dominio completo, n0=12 | 9.965.924 | 0 | 474,47 s |
| Dominio completo, n0=10 | 9.966.320 | 0 | 496,81 s |
| Triangolo base, n0=64 | 2.356 | 0 | <1 s |
| Triangolo base, n0=48 | 1.618 | 0 | <1 s |

Sono passati anche:

- la suite originale di contenimento: 270.000 campioni;
- la suite nuova di copertura esatta e aritmetica a 80 cifre;
- 200 confronti delle identita del campo e dell'Hessiana;
- il caso equilatero con quattro equilibri, un minimo planare e tre selle;
- 2.000 confronti delle formule con una derivazione tramite sistemi lineari;
- 3.015 argomenti e 1.000 intervalli trigonometrici confrontati a 80 cifre;
- 1.260 controlli ad alta precisione degli invarianti e delle derivate automatiche;
- due controlli dell'uscita non riuscita: 12 e 10 box deliberatamente irrisolte,
  entrambi con exit code 1.

Questi test numerici sono controlli diagnostici. La verifica del dominio continuo
e effettuata dalle esecuzioni a intervalli, non dal campionamento.

I risultati riproducibili sono in
[`../certificates/review-2026-10-30/reproduction.json`](../certificates/review-2026-10-30/reproduction.json)
e nei relativi log. I tre JSON nella radice del repository sono copie dei nuovi
risultati. Le esecuzioni precedenti restano nella storia Git e nell'audit locale.

```powershell
python -m pip install -r requirements-review.txt
python reproduce.py --workers 16
```

## Limiti e prossimo controllo esterno

Le due griglie utilizzano la stessa implementazione binary64. Non sono due prove
aritmetiche indipendenti. I confronti mpmath sono test puntuali, non una seconda
certificazione del dominio. I file JSON sono riepiloghi delle esecuzioni:
non contengono l'albero completo delle box per un verificatore separato.

Il passo con maggior valore scientifico dopo questa revisione e una verifica
esterna con due obiettivi distinti:

1. Ricostruire gli argomenti geometrici delle Sezioni 2-5: grado di Z sul dominio
   con vertici risolti, winding, esclusione di due doppi zeri, cuspidi, rotazione
   globale e passaggio a tutti i parametri delle cariche.
2. Ricertificare le formule con una seconda implementazione rigorosa, per esempio
   Arb/MPFI, includendo copertura, simmetrie e regimi singolari; in alternativa,
   progettare un formato delle foglie e un verificatore indipendente piu piccolo.

E stato creato uno snapshot locale di sorgenti, dipendenze dichiarate, log e
manoscritto, con manifest SHA256, nell'audit `2026-10-30-followup`:
`maxwell3-reviewed-2026-10-31.zip`. Prima della diffusione pubblica va associato
a una revisione pubblica identificabile. Le modifiche di questa sessione non
sono state pubblicate. Il risultato puo essere presentato come proposta di prova
assistita dal calcolatore, dichiarando esplicitamente lo stato della revisione.

## File del manoscritto

La versione aggiornata e `paper.pdf`, generata da `paper.html`.
Istruzioni di build in [README.md](README.md).

*Nota successiva (31 ottobre 2026).* Durante la revisione il manoscritto era stato
salvato come `paper-reviewed.pdf`, perche `paper.pdf` era bloccato da un'applicazione.
Integrando la revisione nel repository, le due versioni sono state riunite in `paper.pdf`.
Sono state inoltre precisate la citazione di Martinez-Maure e la dichiarazione sull'uso
dell'IA. L'esecuzione n0=12 e stata ripetuta in una sessione separata con conteggi identici;
log e JSON sono in `../certificates/rerun-2026-10-31/`. Il modello usato per questa revisione
in OpenAI Codex era GPT-6-Astra, secondo quanto riferito dall'autore.

*Nota successiva (10 novembre 2026).* Una seconda verifica, indipendente da questa, e stata
eseguita il 9 novembre 2026 da Claude Fable 5.2: il rapporto e in
[REVIEW-2026-11-09.md](REVIEW-2026-11-09.md). Ha portato a due correzioni del codice
(controllo del dominio di seno e coseno; rimozione di un vecchio percorso non rigoroso) e a
una seconda implementazione, ora in `../independent/`. Le certificazioni sono state
rieseguite con i sorgenti corretti: i conteggi delle box sono identici a quelli della
tabella qui sopra, i tempi sono diversi. Le esecuzioni correnti sono in
`../certificates/review-2026-11-10/`; quelle descritte in questo rapporto restano in
`../certificates/review-2026-10-30/` e si riferiscono ai sorgenti del commit `1d7c1f8`.
