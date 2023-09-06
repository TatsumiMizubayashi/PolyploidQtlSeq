namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// 入力VCFファイル
    /// </summary>
    internal class InputVcf
    {
        /// <summary>
        /// 入力ファイルを作成する。
        /// </summary>
        /// <param name="filePath">VCFファイルPath</param>
        public InputVcf(string filePath)
        {
            if (string.IsNullOrWhiteSpace(filePath)) throw new ArgumentNullException(nameof(filePath));

            Path = System.IO.Path.GetFullPath(filePath);
        }

        /// <summary>
        /// VCFファイルパスを取得する。
        /// </summary>
        internal string Path { get; }
    }
}
