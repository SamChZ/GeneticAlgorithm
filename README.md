# NaCo Opdracht 21/22 – Genetisch Algoritme & Cellulaire Automaat

Gemaakt voor de NaCo opdracht van studiejaar 2021/2022 aan de Universiteit Leiden, groep 16.

Het project heeft twee onderdelen: een genetisch algoritme (GA) en een cellulaire automaat (CA). Het idee is dat het GA het inverse probleem van de CA oplost — je geeft een eindtoestand Cᵗ mee en het GA probeert de bijbehorende begintoestand C⁰ te vinden.

## Uitvoeren

```bash
pip install -r requirements.txt
python implementation.py
```

Wil je alleen het GA testen op een los probleem:

```bash
python test_algorithm.py implementation.py
```

## Wat er in zit

Het GA heeft een paar varianten voor elk onderdeel die je kunt combineren:

- **Selectie:** rank of tournament
- **Crossover:** uniform of one-point
- **Mutatie:** bitstring, bitflip of swap
- **CA:** werkt met binaire (k=2) en ternaire (k=3) representaties
- **Similariteit:** Hamming distance en Levenshtein distance

Bij het uitvoeren van `main()` worden alle combinaties automatisch doorgelopen over 10 problemen, elk 5 keer.

## Benodigde pakketten

- ioh
- NumPy
- pandas
- python-Levenshtein
