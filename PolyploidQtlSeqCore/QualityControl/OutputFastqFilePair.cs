using PolyploidQtlSeqCore.IO;

namespace PolyploidQtlSeqCore.QualityControl
{
    /// <summary>
    /// 出力Fastqファイルペア
    /// </summary>
    internal class OutputFastqFilePair
    {
        /// <summary>
        /// 出力Fastqファイルペアを作成する。
        /// </summary>
        /// <param name="outDirPath">出力ディレクトリPATH</param>
        /// <param name="fastqFilePair">Fastqファイルペア</param>
        public OutputFastqFilePair(string outDirPath, FastqFilePair fastqFilePair)
        {
            if (string.IsNullOrWhiteSpace(outDirPath)) throw new ArgumentException(null, nameof(outDirPath));

            var fastq1Name = Path.GetFileName(fastqFilePair.Fastq1Path);
            Fastq1Path = Path.Combine(outDirPath, fastq1Name);

            var fastq2Name = Path.GetFileName(fastqFilePair.Fastq2Path);
            Fastq2Path = Path.Combine(outDirPath, fastq2Name);

            BaseName = fastqFilePair.BaseName;
        }

        /// <summary>
        /// 出力ファイル1のPATHを取得する。
        /// </summary>
        public string Fastq1Path { get; }

        /// <summary>
        /// 出力ファイル2のPATHを取得する。
        /// </summary>
        public string Fastq2Path { get; }

        /// <summary>
        /// Fastqファイルのベース名を取得する。
        /// </summary>
        public string BaseName { get; }

        /// <summary>
        /// Fastpの引数に変換する。
        /// </summary>
        /// <returns>Fastp引数</returns>
        public string ToFastpArg()
        {
            return $"-o {Fastq1Path} -O {Fastq2Path}";
        }
    }
}
