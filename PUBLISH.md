# Publishing checklist

Four commands put the code on GitHub; the web interface then goes live as a
Hugging Face Space built from the same repository.

## 0. Decide the copyright holder

`LICENSE` and `CITATION.cff` currently name **Belle Damodara Shenoy** under the
MIT licence. If the work belongs to MycoAsia or CSIR-NIO instead, change it now
rather than after the first public commit:

```bash
python deploy/configure.py --github-user belleshenoy1978 --copyright "CSIR-NIO"
```

Leave `--copyright` off to keep the current holder. The same command fills in the
GitHub URLs everywhere they appear, so run it either way.

## 1. GitHub

```bash
cd dna-barcoding
python deploy/configure.py --github-user belleshenoy1978
git init -b main
git add -A
git commit -m "DNA-barcode identification and phylogeny pipeline"
git remote add origin https://github.com/belleshenoy1978/dna-barcoding.git
git push -u origin main
```

Create the empty repository on <https://github.com/new> first (no README, no
.gitignore - this repo already has both). Tests run automatically on Linux and
Windows after the first push; the badge-free workflow is in
`.github/workflows/ci.yml`.

## 2. Hugging Face Space

Full detail in [`deploy/huggingface/DEPLOY.md`](deploy/huggingface/DEPLOY.md).
In short: create a **Docker** Space at <https://huggingface.co/new-space>, then

```bash
git clone https://huggingface.co/spaces/<owner>/barcode-pipeline hf-space
cd hf-space
cp -r ../dna-barcoding/* .
cp deploy/huggingface/README.md README.md     # carries the Space front-matter
git add -A && git commit -m "Barcode pipeline web interface" && git push
```

The first build takes a few minutes. The Space starts in the restricted public
profile: identification disabled, uploads capped, jobs deleted after six hours.

## 3. Your own subdomain (optional)

At your DNS provider, add one record:

| Type | Name | Value |
|---|---|---|
| CNAME | `barcode` | `<owner>-barcode-pipeline.hf.space` |

Then add `barcode.belleshenoy.in` under the Space's **Settings - Custom domain**.
Your main site at `www.belleshenoy.in` is unaffected; the subdomain resolves
independently.

## 4. Afterwards

* Add a screenshot of the running interface to `README.md` - the one asset that
  could not be produced without a browser.
* Tag a release (`git tag v1.0.0 && git push --tags`) so `CITATION.cff` points at
  something citable.
* If you later want identification enabled on the public instance, set
  `BARCODE_ALLOW_BLAST=1` in the Space variables - but read the note in
  `DEPLOY.md` about NCBI rate-limiting first.
