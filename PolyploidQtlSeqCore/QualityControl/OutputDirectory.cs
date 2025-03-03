using PolyploidQtlSeqCore.IO;

namespace PolyploidQtlSeqCore.QualityControl
{
    /// <summary>
    /// 出力ディレクトリ
    /// </summary>
    internal class OutputDirectory
    {
        /// <summary>
        /// 出力ディレクトリ インスタンスを作成する。
        /// </summary>
        /// <param name="outputDirPath">出力ディレクトリPath</param>
        public OutputDirectory(string outputDirPath)
        {
            ArgumentException.ThrowIfNullOrEmpty(outputDirPath);

            Path = System.IO.Path.GetFullPath(outputDirPath);
        }

        /// <summary>
        /// 出力ディレクトリPATHを取得する。
        /// </summary>
        internal string Path { get; }

        /// <summary>
        /// 出力ディレクトリ内のファイルPathを作成する。
        /// </summary>
        /// <param name="fileNames">ファイル名</param>
        /// <returns>出力ディレクトリ内ファイルPATH</returns>
        public string CreateFilePath(params string[] fileNames)
        {
            var items = new[] { Path }.Concat(fileNames).ToArray();
            return System.IO.Path.Combine(items);
        }

        /// <summary>
        /// サブディレクトリを作成し、そのPathを返す。
        /// </summary>
        /// <param name="subDirName">サブディレクトリ名</param>
        /// <returns>サブディレクトリのPATH</returns>
        public string CreateSubDir(string subDirName)
        {
            Create();

            var subDirPath = System.IO.Path.Combine(Path, subDirName);
            if (Directory.Exists(subDirPath)) return subDirPath;

            Directory.CreateDirectory(subDirPath);

            return subDirPath;
        }

        /// <summary>
        /// 出力ディレクトリを作成する。
        /// </summary>
        internal void Create()
        {
            if (Directory.Exists(Path)) return;

            Directory.CreateDirectory(Path);
        }


        /// <summary>
        /// 出力Fastqファイルペアに変換する。
        /// </summary>
        /// <param name="fastqFilePair">Fastqファイルペア</param>
        /// <returns>出力ファイルペア</returns>
        internal OutputFastqFilePair ToOutputFastqFilePair(FastqFilePair fastqFilePair)
        {
            return new OutputFastqFilePair(Path, fastqFilePair);
        }

    }
}
