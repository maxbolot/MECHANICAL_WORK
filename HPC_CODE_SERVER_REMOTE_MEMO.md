# HPC Remote Development Memo: iPad + ThinkPad + SLURM Compute Node

This guide details the procedure to run **`code-server`** on a dynamically allocated **SLURM compute node** on an HPC cluster, tunnel it through your **ThinkPad** over **Tailscale**, and work seamlessly from an **iPad** (with full LaTeX Workshop preview, language tools, and Antigravity AI pair programming).

---

## 1. Architecture Overview

```
┌────────────────────────────────────────────────────────┐
│                      iPad (Client)                     │
│  - Safari Tab 1: code-server (Editor, Git, LaTeX PDF)  │
│  - Safari Tab 2: antigravity.google.com (AI Agent)    │
└──────────────────────────┬─────────────────────────────┘
                           │ Tailscale (Encrypted Mesh)
                           ▼
┌────────────────────────────────────────────────────────┐
│                   ThinkPad P51 (Bridge)                │
│  - Tailscale serve / port forwarder                    │
│  - SSH tunnel to HPC Login Node                        │
│  - Optional: agy daemon / sshfs mount                  │
└──────────────────────────┬─────────────────────────────┘
                           │ SSH Port Forward (-L)
                           ▼
┌────────────────────────────────────────────────────────┐
│                     HPC Login Node                     │
│  - Resolves internal compute node addresses           │
└──────────────────────────┬─────────────────────────────┘
                           │ High-speed internal cluster fabric
                           ▼
┌────────────────────────────────────────────────────────┐
│             SLURM Compute Node (Execution)             │
│  - code-server running on allocated CPUs/RAM           │
│  - Project files, build tools, SLURM jobs              │
│  - Optional: agy --remote-control (if proxy available) │
└────────────────────────────────────────────────────────┘
```

---

## 2. Step 1: Launch `code-server` on the Compute Node

Compute nodes must be allocated through SLURM. You can run interactively with `salloc` or submit a standing batch job with `sbatch`.

### Option A: Interactive Session (`salloc`)
SSH to the HPC login node from the ThinkPad:
```bash
# 1. Allocate resources (e.g. 1 node, 4 CPUs, 16GB RAM for 8 hours)
salloc --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=16G --time=08:00:00

# 2. Check the compute node hostname
NODE=$(hostname)
echo "Running on compute node: $NODE"

# 3. Start code-server on an internal port (e.g. 8888)
code-server --bind-addr 0.0.0.0:8888 --auth none
```

---

### Option B: Batch Job Script (`sbatch`) — Recommended
Create a script `~/start_code_server.slurm` in your HPC home directory:

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --job-name=code-server
#SBATCH --output=%x-%j.out

PORT=8888
NODE=$(hostname)

# Save connection details to a known file
echo "${NODE}:${PORT}" > ~/.code-server-active-node

echo "=========================================================="
echo "code-server started on node: $NODE at port: $PORT"
echo "Submitted SLURM Job ID: $SLURM_JOB_ID"
echo ""
echo "To tunnel from ThinkPad, run:"
echo "ssh -N -L 8080:${NODE}:${PORT} <username>@<hpc-login-host>"
echo "=========================================================="

code-server --bind-addr 0.0.0.0:${PORT} --auth none
```

Submit it:
```bash
sbatch ~/start_code_server.slurm
```
Check status and node assignment:
```bash
cat ~/.code-server-active-node
# Example output: della-i14g2:8888
```

---

## 3. Step 2: Establish the SSH Tunnel on the ThinkPad

SSH can resolve and route internal hostnames through the login node.

On the **ThinkPad terminal**:
```bash
# Forward local ThinkPad port 8080 to compute node port 8888 via login node
ssh -N -L 8080:<COMPUTE_NODE_HOSTNAME>:8888 <username>@<hpc-login-host>
```
*(Example: `ssh -N -L 8080:della-i14g2:8888 max@della.princeton.edu`)*

### Quick Helper Script on ThinkPad
You can automate this with a small bash script `~/bin/connect-hpc-vscode.sh`:
```bash
#!/bin/bash
HPC_USER="max"
HPC_LOGIN="della.princeton.edu"

# Fetch active node from HPC
REMOTE_TARGET=$(ssh ${HPC_USER}@${HPC_LOGIN} "cat ~/.code-server-active-node 2>/dev/null")

if [ -z "$REMOTE_TARGET" ]; then
  echo "Error: No active code-server found in ~/.code-server-active-node on HPC."
  exit 1
fi

echo "Connecting to $REMOTE_TARGET via $HPC_LOGIN..."
ssh -N -L 8080:${REMOTE_TARGET} ${HPC_USER}@${HPC_LOGIN}
```

---

## 4. Step 3: Expose via Tailscale & Connect from iPad

### Expose on ThinkPad
Make the local forwarded port accessible over your Tailscale tailnet:
```bash
tailscale serve --bg http://127.0.0.1:8080
```
*(Or for HTTPS with valid TLS certificates: `tailscale serve https / http://127.0.0.1:8080`)*

### Open on iPad
Open Safari on your iPad and go to:
```text
http://<your-thinkpad-name>.tailnet-name.ts.net:8080
# Example:
http://max-thinkpad-p51.tail912b28.ts.net:8080
```

> [!TIP]
> **iPad Home Screen App**: In Safari on iPad, tap the **Share** button -> **Add to Home Screen**. This runs `code-server` in standalone fullscreen mode without the browser URL bar or navigation tabs.

---

## 5. Step 4: Antigravity AI (`agy`) Integration

You can integrate Antigravity pair programming into this workflow in two ways:

### Method A: Direct on Compute Node (If HPC has Outbound Internet / Proxy)
Most clusters provide outbound HTTP access via an environment proxy:
1. Inside `code-server` on your iPad, open an integrated terminal:
   ```bash
   # If required by your cluster, export HTTP proxies
   export http_proxy=http://proxy.cluster:8080
   export https_proxy=http://proxy.cluster:8080

   # Launch agy in remote control mode
   agy --remote-control
   ```
2. Open **`antigravity.google.com`** in a second Safari tab on your iPad.
3. You can prompt the AI agent from the web interface, and it edits files and compiles code directly inside the HPC compute node.

---

### Method B: Hybrid SSHFS (If Compute Nodes Have No Outbound Internet)
If the HPC compute nodes are completely air-gapped from the internet:
1. On the **ThinkPad**, mount the HPC project directory via `sshfs`:
   ```bash
   mkdir -p ~/mnt/hpc_project
   sshfs <username>@<hpc-login-host>:/path/to/project ~/mnt/hpc_project \
     -o reconnect,ServerAliveInterval=15,ServerAliveCountMax=3
   ```
2. Run `agy --remote-control` on the **ThinkPad** inside `~/mnt/hpc_project`.
3. The AI agent runs on the ThinkPad with full Google Cloud internet access, modifying files over the `sshfs` mount.
4. Any change made by the agent is immediately reflected in `code-server` on the compute node and on your iPad screen.

---

## 6. Quick Cheat Sheet

| Action | Machine | Command |
|---|---|---|
| **Start Job** | HPC Login | `sbatch ~/start_code_server.slurm` |
| **Check Node** | HPC Login | `cat ~/.code-server-active-node` |
| **Tunnel Port** | ThinkPad | `ssh -N -L 8080:<node>:8888 user@hpc-login` |
| **Tailscale Serve** | ThinkPad | `tailscale serve --bg http://127.0.0.1:8080` |
| **Open UI** | iPad Safari | `http://max-thinkpad-p51...:8080` |
| **AI Agent** | iPad / Browser | `https://antigravity.google.com` |
