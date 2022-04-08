all:
	edictor wordlist --dataset=cldf-datasets/carvalhopurus/cldf/cldf-metadata.json --preprocessing=pkg/preprocessing.py --name=data/carvalhopurus
	edictor wordlist --dataset=cldf-datasets/meloniromance/cldf/cldf-metadata.json --addon=cognacy:cogids --name=data/meloniromance
	edictor wordlist --dataset=cldf-datasets/wangbai/cldf/cldf-metadata.json --preprocessing=pkg/preprocessing.py --name=data/wangbai
	edictor wordlist --dataset=cldf-datasets/hillburmish/cldf/cldf-metadata.json --addon=partial_cognacy:cogids --name=data/hillburmish
	edictor wordlist --dataset=cldf-datasets/yanglalo/cldf/cldf-metadata.json --preprocessing=pkg/preprocessing.py --name=data/yanglalo
	edictor wordlist --dataset=cldf-datasets/ltkkaren/cldf/cldf-metadata.json --addon=partial_cognacy:cogids --name=data/ltkkaren
