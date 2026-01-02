import React, { useState } from 'react';

interface Atom {
  id: number;
  element: string;
  position: [number, number, number];
  label: string;
  charge: number;
}

interface Bond {
  from: number;
  to: number;
}

interface ModelInfoProps {
  fileName: string;
  atoms: Atom[];
  bonds: Bond[];
  onSelectAtom?: (atomId: number) => void;
}

export default function ModelInfo({ fileName, atoms, bonds, onSelectAtom }: ModelInfoProps): React.ReactElement {
  const [activeTab, setActiveTab] = useState<'summary' | 'atoms' | 'bonds'>('summary');
  
  // Calculate some useful statistics
  const elementCounts = atoms.reduce((counts, atom) => {
    counts[atom.element] = (counts[atom.element] || 0) + 1;
    return counts;
  }, {} as Record<string, number>);
  
  const sortedElements = Object.entries(elementCounts).sort((a, b) => b[1] - a[1]);
  
  const averageCharge = atoms.length > 0
    ? atoms.reduce((sum, atom) => sum + atom.charge, 0) / atoms.length
    : 0;
    
  const totalCharge = atoms.reduce((sum, atom) => sum + atom.charge, 0);
  
  // Calculate min and max coordinates for each axis
  let minX = Infinity, maxX = -Infinity;
  let minY = Infinity, maxY = -Infinity;
  let minZ = Infinity, maxZ = -Infinity;
  
  atoms.forEach(atom => {
    minX = Math.min(minX, atom.position[0]);
    maxX = Math.max(maxX, atom.position[0]);
    minY = Math.min(minY, atom.position[1]);
    maxY = Math.max(maxY, atom.position[1]);
    minZ = Math.min(minZ, atom.position[2]);
    maxZ = Math.max(maxZ, atom.position[2]);
  });
  
  const dimensions = {
    x: maxX - minX,
    y: maxY - minY,
    z: maxZ - minZ,
  };

  return (
    <div className="model-info-panel">
      <div className="model-file-info">
        <h2>MODEL: {fileName}</h2>
      </div>
      
      <div className="tabs">
        <button 
          className={activeTab === 'summary' ? 'active' : ''} 
          onClick={() => setActiveTab('summary')}
        >
          Summary
        </button>
        <button 
          className={activeTab === 'atoms' ? 'active' : ''} 
          onClick={() => setActiveTab('atoms')}
        >
          Atoms ({atoms.length})
        </button>
        <button 
          className={activeTab === 'bonds' ? 'active' : ''} 
          onClick={() => setActiveTab('bonds')}
        >
          Bonds ({bonds.length})
        </button>
      </div>
      
      <div className="tab-content">
        {activeTab === 'summary' && (
          <div className="summary-tab">
            <div className="stats-group">
              <h3>Composition</h3>
              <div className="element-list">
                {sortedElements.map(([element, count]) => (
                  <div key={element} className="element-item">
                    <span className="element-symbol">{element}</span>
                    <span className="element-count">{count}</span>
                  </div>
                ))}
              </div>
            </div>
            
            <div className="stats-group">
              <h3>Dimensions (Å)</h3>
              <table className="stats-table">
                <tbody>
                  <tr>
                    <td>X-axis:</td>
                    <td>{dimensions.x.toFixed(2)}</td>
                  </tr>
                  <tr>
                    <td>Y-axis:</td>
                    <td>{dimensions.y.toFixed(2)}</td>
                  </tr>
                  <tr>
                    <td>Z-axis:</td>
                    <td>{dimensions.z.toFixed(2)}</td>
                  </tr>
                </tbody>
              </table>
            </div>
            
            <div className="stats-group">
              <h3>Charges</h3>
              <table className="stats-table">
                <tbody>
                  <tr>
                    <td>Average Charge:</td>
                    <td>{averageCharge.toFixed(4)}</td>
                  </tr>
                  <tr>
                    <td>Total Charge:</td>
                    <td>{totalCharge.toFixed(4)}</td>
                  </tr>
                </tbody>
              </table>
            </div>
          </div>
        )}
        
        {activeTab === 'atoms' && (
          <div className="atoms-tab">
            <div className="atom-list-container">
              <table className="atom-table">
                <thead>
                  <tr>
                    <th>ID</th>
                    <th>Element</th>
                    <th>Label</th>
                    <th>Position (X, Y, Z)</th>
                    <th>Charge</th>
                  </tr>
                </thead>
                <tbody>
                  {atoms.map((atom) => (
                    <tr 
                      key={atom.id} 
                      onClick={() => onSelectAtom && onSelectAtom(atom.id)}
                      className="atom-row"
                    >
                      <td>{atom.id}</td>
                      <td>{atom.element}</td>
                      <td>{atom.label}</td>
                      <td>
                        ({atom.position[0].toFixed(3)}, {atom.position[1].toFixed(3)}, {atom.position[2].toFixed(3)})
                      </td>
                      <td>{atom.charge.toFixed(3)}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        )}
        
        {activeTab === 'bonds' && (
          <div className="bonds-tab">
            <div className="bond-list-container">
              <table className="bond-table">
                <thead>
                  <tr>
                    <th>Bond</th>
                    <th>From Atom</th>
                    <th>To Atom</th>
                  </tr>
                </thead>
                <tbody>
                  {bonds.map((bond, index) => (
                    <tr key={index} className="bond-row">
                      <td>{index + 1}</td>
                      <td>{bond.from}</td>
                      <td>{bond.to}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        )}
      </div>
      
      <style>{`
        .model-info-panel {
          width: 100%;
          height: 100%;
          padding: 8px 10px;
          color: #7AFFB2;
          overflow: hidden;
          display: flex;
          flex-direction: column;
          font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
          font-size: 11px;
        }

        .model-file-info {
          margin-bottom: 8px;
        }

        .model-file-info h2 {
          margin: 0 0 2px;
          font-size: 12px;
          color: #07DBFA;
          letter-spacing: 1px;
          text-shadow: 0 0 5px rgba(7, 219, 250, 0.5);
        }

        .model-file-info p {
          margin: 0;
          font-size: 10px;
          color: #7AFFB2;
          opacity: 0.8;
        }

        .tabs {
          display: flex;
          border-bottom: 1px solid rgba(7, 219, 250, 0.3);
          margin-bottom: 8px;
        }

        .tabs button {
          padding: 4px 8px;
          background: none;
          border: none;
          border-bottom: 2px solid transparent;
          cursor: pointer;
          font-size: 10px;
          font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
          color: #A6A8F8;
          transition: all 0.2s;
          text-transform: uppercase;
          letter-spacing: 0.5px;
        }

        .tabs button:hover {
          color: #07DBFA;
        }

        .tabs button.active {
          color: #07DBFA;
          border-bottom: 2px solid #07DBFA;
        }

        .tab-content {
          padding: 4px 0;
          flex: 1;
          overflow-y: auto;
          max-height: 280px;
        }

        .stats-group {
          margin-bottom: 10px;
        }

        .stats-group h3 {
          margin: 0 0 4px;
          font-size: 11px;
          color: #F8D2A6;
          border-bottom: 1px dotted rgba(248, 210, 166, 0.3);
          padding-bottom: 2px;
          text-transform: uppercase;
        }

        .element-list {
          display: flex;
          flex-wrap: wrap;
          gap: 6px;
        }

        .element-item {
          display: flex;
          align-items: center;
          gap: 4px;
          padding: 2px 6px;
          background-color: rgba(7, 219, 250, 0.1);
          border-radius: 3px;
          font-size: 11px;
          border: 1px solid rgba(7, 219, 250, 0.3);
        }

        .element-symbol {
          font-weight: bold;
        }

        .stats-table {
          width: 100%;
          border-collapse: collapse;
          font-size: 11px;
        }

        .stats-table td {
          padding: 3px 8px;
          border-bottom: 1px solid rgba(7, 219, 250, 0.2);
        }

        .stats-table td:first-child {
          font-weight: 500;
          color: #F8D2A6;
        }

        .atom-list-container, .bond-list-container {
          max-height: 280px;
          overflow-y: auto;
        }

        .atom-table, .bond-table {
          width: 100%;
          border-collapse: collapse;
          font-size: 10px;
        }

        .atom-table th, .bond-table th {
          background-color: rgba(0, 0, 0, 0.4);
          padding: 4px 6px;
          text-align: left;
          position: sticky;
          top: 0;
          color: #07DBFA;
          font-weight: normal;
          text-transform: uppercase;
          letter-spacing: 0.5px;
          font-size: 10px;
        }

        .atom-row, .bond-row {
          cursor: pointer;
          transition: background-color 0.2s;
        }

        .atom-row:hover, .bond-row:hover {
          background-color: rgba(122, 255, 178, 0.1);
        }

        .atom-table td, .bond-table td {
          padding: 2px 6px;
          border-bottom: 1px solid rgba(7, 219, 250, 0.2);
        }
      `}</style>
    </div>
  );
}