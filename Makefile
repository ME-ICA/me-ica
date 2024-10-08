
build_dev:
	docker build -t devcontainer .

run_dev:
	docker run --name devcontainer -dt devcontainer 

exec_dev:
	docker exec -it devcontainer /bin/bash
