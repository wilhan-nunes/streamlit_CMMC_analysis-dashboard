server-compose-build-nocache:
	docker-compose --compatibility build --no-cache

server-compose-build:
	docker-compose --compatibility --verbose build 

server-compose-interactive:
	docker-compose --compatibility build
	docker-compose --compatibility -f docker-compose.yml -f docker-compose-dev.yml up

server-compose:
	docker-compose --compatibility build
	docker-compose --compatibility -f docker-compose.yml -f docker-compose-dev.yml up -d

server-compose-production:
	docker-compose --compatibility build
	docker-compose --compatibility -f docker-compose.yml up -d

attach:
	docker exec -i -t streamlit-cmmc-dashboard /bin/bash

init_modules:
	git submodule update --init --recursive