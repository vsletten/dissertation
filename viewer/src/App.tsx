import { useState } from 'react';
import FileUploader from './components/FileUploader';
import MoleculeViewer from './components/MoleculeViewer';
import ModelInfo from './components/ModelInfo';
import { parseModelFile } from './utils/cerius2Parser';

function App() {
  const [fileName, setFileName] = useState<string>('');
  const [atoms, setAtoms] = useState<Array<{
    id: number;
    element: string;
    position: [number, number, number];
    label: string;
    charge: number;
  }>>([]);
  const [bonds, setBonds] = useState<Array<{ from: number; to: number }>>([]);
  const [error, setError] = useState<string | null>(null);

  const handleFileLoaded = (content: string, name: string) => {
    try {
      setError(null);
      
      // Parse the model file to get atoms and bonds directly
      const result = parseModelFile(content);
      
      if (!result) {
        throw new Error('Failed to parse model file');
      }
      
      setFileName(name);
      setAtoms(result.atoms);
      setBonds(result.bonds);
    } catch (err) {
      setError(`Error processing file: ${err instanceof Error ? err.message : String(err)}`);
      setFileName('');
      setAtoms([]);
      setBonds([]);
    }
  };

  return (
    <div className="app">
      {atoms.length === 0 ? (
        <div className="upload-section">
          <FileUploader onFileLoaded={handleFileLoaded} />
          {error && <div className="error-message">{error}</div>}
          
          <div className="instruction-panel">
            <h2>Getting Started</h2>
            <p style={{ color: '#7AFFB2', fontSize: '0.8rem', marginTop: '5px' }}>
              Upload a Cerius2 model file (.msi, .car) to visualize molecular structures.
            </p>
          </div>
        </div>
      ) : (
        <div className="viewer-section">
          <div className="viewer-container">
            <MoleculeViewer atoms={atoms} bonds={bonds} />
            
            <div className="model-info-overlay">
              <ModelInfo 
                fileName={fileName}
                atoms={atoms}
                bonds={bonds}
              />
              
              <button
                className="new-file-button"
                onClick={() => {
                  setFileName('');
                  setAtoms([]);
                  setBonds([]);
                  setError(null);
                }}
              >
                Load New File
              </button>
            </div>
          </div>
        </div>
      )}
      
      <style>{`
        .app {
          width: 100vw;
          height: 100vh;
          overflow: hidden;
          background-color: #111;
          position: relative;
        }
        
        .upload-section {
          width: 100%;
          height: 100%;
          display: flex;
          flex-direction: column;
          justify-content: center;
          align-items: center;
          padding: 2rem;
        }
        
        .error-message {
          color: #d32f2f;
          background-color: rgba(255, 235, 238, 0.8);
          padding: 10px 15px;
          border-radius: 4px;
          margin-top: 10px;
          font-size: 0.9rem;
        }
        
        .instruction-panel {
          background-color: rgba(10, 20, 30, 0.7);
          border-radius: 4px;
          padding: 15px;
          box-shadow: 0 0 20px rgba(7, 219, 250, 0.2);
          margin-top: 20px;
          max-width: 600px;
          border: 1px solid rgba(7, 219, 250, 0.3);
          font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
        }
        
        .instruction-panel h2 {
          margin-top: 0;
          color: #07DBFA;
          font-size: 1.1rem;
          text-transform: uppercase;
          letter-spacing: 1px;
          text-shadow: 0 0 5px rgba(7, 219, 250, 0.5);
        }
        
        .viewer-section {
          width: 100%;
          height: 100%;
          position: relative;
        }
        
        .viewer-container {
          width: 100%;
          height: 100%;
          position: relative;
        }
        
        .model-info-overlay {
          position: absolute;
          bottom: 20px;
          right: 20px;
          width: 220px;
          max-height: 280px;
          background-color: rgba(10, 20, 30, 0.7);
          backdrop-filter: blur(5px);
          border: 1px solid rgba(7, 219, 250, 0.3);
          border-radius: 4px;
          color: #7AFFB2;
          overflow: auto;
          padding: 0;
          box-shadow: 0 0 15px rgba(7, 219, 250, 0.2),
                     inset 0 0 5px rgba(7, 219, 250, 0.05);
          z-index: 10;
        }
        
        .new-file-button {
          background-color: rgba(7, 219, 250, 0.2);
          color: #07DBFA;
          border: 1px solid rgba(7, 219, 250, 0.4);
          padding: 2px 4px;
          border-radius: 2px;
          font-size: 0.45rem;
          font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
          cursor: pointer;
          transition: all 0.2s;
          margin: 2px 6px 4px;
          width: calc(100% - 16px);
          text-transform: uppercase;
          letter-spacing: 1px;
          height: 16px;
        }
        
        .new-file-button:hover {
          background-color: rgba(7, 219, 250, 0.3);
          box-shadow: 0 0 10px rgba(7, 219, 250, 0.4);
        }
      `}</style>
      
      <style>{`
        * {
          box-sizing: border-box;
          margin: 0;
          padding: 0;
        }
        
        body {
          font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
          font-size: 16px;
          line-height: 1.5;
          color: #eee;
          overflow: hidden;
        }
        
        button {
          cursor: pointer;
        }
      `}</style>
    </div>
  );
}

export default App;