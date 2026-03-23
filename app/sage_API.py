from flask import Flask, jsonify, request
import json
from pathlib import Path
import subprocess

app = Flask(__name__)
DATA_ROOT = Path(__file__) / "data"
            
@app.route('/get_curves/<int:p>/<int:n>', methods=['GET'])
def enumerate_order(p, n):
    """Enumerate curves for a specific number field D over F_{p^n}"""
    try:
        D = int(request.args.get('D'))
        
        # Call standalone script in subprocess to avoid Sage segfaults
        result = subprocess.run(
            ['sage', '-python', 'compute_hcp.py', str(p), str(n), str(D)],
            cwd=str(Path(__file__).parent),
            capture_output=True,
            text=True,
            timeout=30
        )
        
        if result.returncode != 0:
            return jsonify({
                'success': False, 
                'error': 'Computation failed',
                'stderr': result.stderr
            }), 500
        
        data = json.loads(result.stdout)
        return jsonify(data)
        
    except Exception as e:
        import traceback
        return jsonify({'success': False, 'error': str(e), 'trace': traceback.format_exc()}), 500

@app.route('/health', methods=['GET'])
def health():
    return jsonify({'status': 'ok'})

if __name__ == '__main__':
    app.run(host='127.0.0.1', port=5001, debug=False)