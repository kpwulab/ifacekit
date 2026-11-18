Running the ifacekit portal as an internal service
===================================================
1. Install into /your_path/ifacekit (full guide: INSTALL.txt; replace /your_path by the folder you use). On Ubuntu:
       sudo apt install python3-venv libcairo2
       sudo mkdir -p /your_path && sudo chown $USER /your_path
       tar xzf ifacekit.tgz -C /your_path          (the archive unpacks to a folder called ifacekit)
       cd /your_path/ifacekit && python3 -m venv .venv && . .venv/bin/activate
       pip install ".[ligand,pdf,server]"
   (libcairo2 is needed for the PNG and PDF files; without it the portal still works but offers SVG only.)
2. Try it:   server/ifacekit-server.sh --host 127.0.0.1 --port 8000      then open http://127.0.0.1:8000
   Check:    curl localhost:8000/healthz     (version, analyses running and waiting, stored jobs, free disk, local Mol*)
3. As a service: server/deploy/ifacekit.service (systemd). Create the user and the data folder first:
       sudo useradd --system --home /var/lib/ifacekit ifacekit ; sudo mkdir -p /var/lib/ifacekit/jobs ; sudo chown -R ifacekit: /var/lib/ifacekit
4. Behind nginx with TLS (optional on a private network): server/deploy/nginx-ifacekit.conf.example
5. Settings (environment variables, see the service file):
       IFACEKIT_DATA           job folder (uploads and results)
       IFACEKIT_KEEP_DAYS      jobs older than this many days are deleted, checked every IFACEKIT_PURGE_MINUTES (default 7 days, 60 minutes)
       IFACEKIT_MAX_JOBS       analyses at the same time (default 4); more requests wait up to IFACEKIT_QUEUE_WAIT seconds (600)
       IFACEKIT_MEM_MB         memory limit per analysis process in MB (default 8192, 0 = none); uses ulimit -v, Linux
       IFACEKIT_THREADS        server threads (default 16)
   Keep one gunicorn worker process (the launcher does): the limit on simultaneous analyses is counted per process.
6. Mol* for the 3D pages is served from server/static/molstar (molstar.js, molstar.css, LICENSE of Mol*, MIT), so the 3D pages work without internet.
   The portal still needs outgoing HTTPS to files.rcsb.org, www.ebi.ac.uk and pdbj.org for "fetch by PDB ID".
7. Updates: stop the service, unpack the new archive over the folder (keep .venv), pip install -U ., restart.
8. Look after the service:  sudo systemctl status ifacekit   |   journalctl -u ifacekit -f   |   sudo systemctl restart ifacekit
