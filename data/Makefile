download:
	mkdir ldgm
	curl -L -o - https://www.dropbox.com/s/km9qe0q4lswsq25/AFR.tar.gz | tar -xz -C ldgm
	curl -L -o - https://www.dropbox.com/s/3g8x1rkj7mqnypr/AMR.tar.gz | tar -xz -C ldgm
	curl -L -o - https://www.dropbox.com/s/v19rm6s9or47zgr/EAS.tar.gz | tar -xz -C ldgm
	curl -L -o - https://www.dropbox.com/s/p5l2l2peny9ox3n/EUR.tar.gz | tar -xz -C ldgm
	curl -L -o - https://www.dropbox.com/s/izv0sf7yunc4uhb/SAS.tar.gz | tar -xz -C ldgm
	curl -L -o - https://www.dropbox.com/s/sk7f3zvtkjluh2m/sumstats.tar.gz | tar -xz
	mv ldgm/*/* ldgm
	ls -d */ | xargs rmdir

download_all:
        curl -L -o ldgm.tar.gz https://www.dropbox.com/sh/1huaxgad2bjjv9a/AACI6pRZpr6gmK5E1mo4fck6a
        tar -xvf ldgm.tar.gz

