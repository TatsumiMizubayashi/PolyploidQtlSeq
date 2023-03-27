using PolyploidQtlSeqCore.IO;

namespace PolyploidQtlSeqCore.QualityControl
{
    /// <summary>
    /// Fastp HTMLレポートファイル
    /// </summary>
    internal class FastpHtmlReportFile
    {
        private const string REPORT_DIR_NAME = "Report";

        /// <summary>
        /// Fastp HTMLレポートファイルを作成する。
        /// </summary>
        /// <param name="outputDir">出力ディレクトリ</param>
        /// <param name="inputFilePair">入力ファイルペア</param>
        public FastpHtmlReportFile(OutputDirectory outputDir, FastqFilePair inputFilePair)
        {
            outputDir.CreateSubDir(REPORT_DIR_NAME);
            Value = outputDir.CreateFilePath(REPORT_DIR_NAME, inputFilePair.BaseName + ".html");
        }

        /// <summary>
        /// HTMLレポートファイルPathを取得する。
        /// </summary>
        public string Value { get; }

        /// <summary>
        /// Fastpの引数に変換する。
        /// </summary>
        /// <returns>Fastp引数</returns>
        public string ToFastpArg()
        {
            return $"-h {Value}";
        }
    }
}
