# Docker Quick Start Guide

## TL;DR - Get Started in 3 Steps

```bash
# 1. Build the image
docker build -t cptac-browser .

# 2. Run the container
docker run -d -p 8000:8000 -v cptac-data:/root/.cptac --name cptac-browser cptac-browser

# 3. Open browser
# Navigate to http://localhost:8000
```

## Using Docker Compose (Even Easier)

```bash
# Start everything
docker-compose up -d

# View logs
docker-compose logs -f

# Stop everything
docker-compose down
```

## Common Docker Commands

### Building

```bash
# Build image
docker build -t cptac-browser .

# Build with no cache (clean build)
docker build --no-cache -t cptac-browser .

# Build with Docker Compose
docker-compose build
```

### Running

```bash
# Run in foreground (see output)
docker run -p 8000:8000 cptac-browser

# Run in background (detached)
docker run -d -p 8000:8000 --name cptac-browser cptac-browser

# Run with persistent cache
docker run -d -p 8000:8000 -v cptac-data:/root/.cptac --name cptac-browser cptac-browser

# Run on different port
docker run -d -p 8080:8000 --name cptac-browser cptac-browser

# Run with more memory
docker run -d -p 8000:8000 -m 4g --name cptac-browser cptac-browser
```

### Managing Containers

```bash
# List running containers
docker ps

# List all containers (including stopped)
docker ps -a

# Stop container
docker stop cptac-browser

# Start stopped container
docker start cptac-browser

# Restart container
docker restart cptac-browser

# Remove container
docker rm cptac-browser

# Remove running container (force)
docker rm -f cptac-browser
```

### Logs and Debugging

```bash
# View logs
docker logs cptac-browser

# Follow logs in real-time
docker logs -f cptac-browser

# View last 100 lines
docker logs --tail 100 cptac-browser

# Run interactively to see errors
docker run -it -p 8000:8000 cptac-browser

# Execute command in running container
docker exec -it cptac-browser /bin/bash

# Inspect container
docker inspect cptac-browser
```

### Volume Management

```bash
# List volumes
docker volume ls

# Inspect volume
docker volume inspect cptac-data

# Remove volume (clears cache!)
docker volume rm cptac-data

# Remove all unused volumes
docker volume prune
```

### Images

```bash
# List images
docker images

# Remove image
docker rmi cptac-browser

# View image layers
docker history cptac-browser

# Check image size
docker images cptac-browser
```

### Cleanup

```bash
# Remove stopped containers
docker container prune

# Remove unused images
docker image prune

# Remove unused volumes
docker volume prune

# Remove everything unused (nuclear option)
docker system prune -a

# With Docker Compose - stop and remove volumes
docker-compose down -v
```

## Persistent Data

The CPTAC library downloads large datasets (~100MB-2GB per cancer type) on first use. To avoid re-downloading:

### Option 1: Named Volume (Recommended)

```bash
docker run -p 8000:8000 -v cptac-data:/root/.cptac cptac-browser
```

**Pros:** Managed by Docker, survives container removal
**Cons:** Data stored in Docker's directory

### Option 2: Bind Mount

```bash
mkdir -p ./cptac_cache
docker run -p 8000:8000 -v $(pwd)/cptac_cache:/root/.cptac cptac-browser
```

**Pros:** Data accessible in your project directory
**Cons:** Path-dependent, less portable

## Troubleshooting

### Port Already in Use

```bash
# Check what's using port 8000
lsof -i :8000

# Use different port
docker run -p 8080:8000 cptac-browser
```

### Container Exits Immediately

```bash
# Check logs
docker logs cptac-browser

# Run interactively
docker run -it -p 8000:8000 cptac-browser
```

### Out of Memory

```bash
# Increase memory limit
docker run -m 4g -p 8000:8000 cptac-browser

# Check Docker memory settings
docker info | grep Memory
```

### Slow Performance

```bash
# Check resource usage
docker stats cptac-browser

# Use persistent volume to cache data
docker run -v cptac-data:/root/.cptac cptac-browser
```

### Permission Issues

```bash
# Run as current user (Linux)
docker run --user $(id -u):$(id -g) -p 8000:8000 cptac-browser
```

### Can't Connect to Container

```bash
# Check if container is running
docker ps

# Check if port is exposed
docker port cptac-browser

# Check container IP
docker inspect -f '{{range.NetworkSettings.Networks}}{{.IPAddress}}{{end}}' cptac-browser

# Test from inside container
docker exec cptac-browser curl http://localhost:8000
```

## Production Deployment

### Behind Reverse Proxy (nginx)

```nginx
server {
    listen 80;
    server_name cptac.example.com;

    location / {
        proxy_pass http://localhost:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;

        # WebSocket support
        proxy_http_version 1.1;
        proxy_set_header Upgrade $http_upgrade;
        proxy_set_header Connection "upgrade";
    }
}
```

### Auto-Restart on Reboot

```bash
docker run -d --restart unless-stopped -p 8000:8000 cptac-browser
```

Or with Docker Compose (already configured in `docker-compose.yml`):

```yaml
restart: unless-stopped
```

## Resources

- [Docker Documentation](https://docs.docker.com/)
- [Docker Compose Documentation](https://docs.docker.com/compose/)
- [Shiny for Python](https://shiny.posit.co/py/)
