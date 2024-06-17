from variable import ntfy_server

def Generate_cURL_cmd(msg: str) -> str:
    if ntfy_server != "":
        return f"curl -d {msg} {ntfy_server}"
    else:
        return ""
