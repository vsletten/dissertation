import { useState, useRef, ChangeEvent } from 'react';

interface FileUploaderProps {
  onFileLoaded: (content: string, fileName: string) => void;
}

export default function FileUploader({ onFileLoaded }: FileUploaderProps): React.ReactElement {
  const [fileName, setFileName] = useState<string>('');
  const [isLoading, setIsLoading] = useState<boolean>(false);
  const [error, setError] = useState<string | null>(null);
  const fileInputRef = useRef<HTMLInputElement>(null);

  const handleFileChange = (event: ChangeEvent<HTMLInputElement>) => {
    const files = event.target.files;
    if (!files || files.length === 0) {
      return;
    }
    
    const file = files[0];
    setFileName(file.name);
    setIsLoading(true);
    setError(null);
    
    const reader = new FileReader();
    
    reader.onload = (e) => {
      try {
        const content = e.target?.result as string;
        onFileLoaded(content, file.name);
        setIsLoading(false);
      } catch (err) {
        setError(`Failed to parse file: ${err instanceof Error ? err.message : String(err)}`);
        setIsLoading(false);
      }
    };
    
    reader.onerror = () => {
      setError('Failed to read file');
      setIsLoading(false);
    };
    
    reader.readAsText(file);
  };

  const handleDragOver = (event: React.DragEvent<HTMLDivElement>) => {
    event.preventDefault();
  };
  
  const handleDrop = (event: React.DragEvent<HTMLDivElement>) => {
    event.preventDefault();
    
    const files = event.dataTransfer.files;
    if (!files || files.length === 0) {
      return;
    }
    
    const file = files[0];
    setFileName(file.name);
    setIsLoading(true);
    setError(null);
    
    const reader = new FileReader();
    
    reader.onload = (e) => {
      try {
        const content = e.target?.result as string;
        onFileLoaded(content, file.name);
        setIsLoading(false);
      } catch (err) {
        setError(`Failed to parse file: ${err instanceof Error ? err.message : String(err)}`);
        setIsLoading(false);
      }
    };
    
    reader.onerror = () => {
      setError('Failed to read file');
      setIsLoading(false);
    };
    
    reader.readAsText(file);
  };

  const handleButtonClick = () => {
    fileInputRef.current?.click();
  };

  return (
    <div className="file-uploader">
      <div 
        className="uploader-dropzone"
        onDragOver={handleDragOver}
        onDrop={handleDrop}
        onClick={handleButtonClick}
      >
        <input 
          type="file" 
          ref={fileInputRef}
          onChange={handleFileChange}
          accept=".msi,.car,.pdb,.xyz"
          style={{ display: 'none' }}
        />
        <div className="uploader-content">
          {isLoading ? (
            <div className="loading-spinner">Loading...</div>
          ) : (
            <>
              <svg 
                xmlns="http://www.w3.org/2000/svg" 
                width="24" 
                height="24" 
                viewBox="0 0 24 24"
                fill="none" 
                stroke="currentColor" 
                strokeWidth="2" 
                strokeLinecap="round" 
                strokeLinejoin="round"
              >
                <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path>
                <polyline points="17 8 12 3 7 8"></polyline>
                <line x1="12" y1="3" x2="12" y2="15"></line>
              </svg>
              <p>
                {fileName 
                  ? `Selected file: ${fileName}`
                  : 'Drag & drop a Cerius2 file here, or click to select'
                }
              </p>
              <p className="file-format-note">
                Supports .msi, .car, .pdb, and .xyz files
              </p>
            </>
          )}
        </div>
      </div>
      {error && <div className="error-message">{error}</div>}
      <style>{`
        .file-uploader {
          margin-bottom: 20px;
          width: 100%;
          max-width: 600px;
        }
        
        .uploader-dropzone {
          border: 2px dashed rgba(7, 219, 250, 0.5);
          border-radius: 6px;
          padding: 40px 20px;
          text-align: center;
          cursor: pointer;
          transition: all 0.3s;
          background-color: rgba(10, 20, 30, 0.5);
          box-shadow: 0 0 30px rgba(7, 219, 250, 0.2);
        }
        
        .uploader-dropzone:hover {
          border-color: rgba(7, 219, 250, 0.8);
          box-shadow: 0 0 50px rgba(7, 219, 250, 0.3);
          transform: scale(1.02);
        }
        
        .uploader-content {
          display: flex;
          flex-direction: column;
          align-items: center;
          justify-content: center;
          color: #7AFFB2;
          font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
        }
        
        .uploader-content svg {
          margin-bottom: 20px;
          color: #07DBFA;
          stroke-width: 1.5;
          filter: drop-shadow(0 0 8px rgba(7, 219, 250, 0.5));
        }
        
        .file-format-note {
          font-size: 0.75rem;
          color: rgba(166, 168, 248, 0.8);
          margin-top: 10px;
          letter-spacing: 0.5px;
        }
        
        .error-message {
          color: #ff4d4d;
          margin-top: 15px;
          font-size: 0.8rem;
          text-shadow: 0 0 10px rgba(255, 77, 77, 0.3);
          font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
        }
        
        .loading-spinner {
          display: flex;
          align-items: center;
          justify-content: center;
          color: #07DBFA;
          font-size: 1rem;
          letter-spacing: 0.5px;
        }
      `}</style>
    </div>
  );
}