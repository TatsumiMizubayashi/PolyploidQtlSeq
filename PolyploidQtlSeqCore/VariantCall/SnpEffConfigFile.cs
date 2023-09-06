namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// SnpEffコンフィグファイル
    /// </summary>
    internal class SnpEffConfigFile
    {
        /// <summary>
        /// SnpEffConfigFileを指定する。
        /// 既定のSnpEffConfigファイルを使用する場合は空文字を指定する。
        /// </summary>
        /// <param name="filePath">snpEffConfigファイルのPath</param>
        public SnpEffConfigFile(string filePath)
        {
            HasFile = !string.IsNullOrEmpty(filePath);

            if(HasFile)
            {
                if (!File.Exists(filePath)) throw new FileNotFoundException(filePath);
                Path = System.IO.Path.GetFullPath(filePath);
            }
            else
            {
                Path = "";
            }
        }

        /// <summary>
        /// SnpEffConfigファイルのPathを指定する。
        /// </summary>
        internal string Path { get; }

        /// <summary>
        /// SnpEffコンフィグファイルが指定されているかどうか。
        /// </summary>
        internal bool HasFile { get; }

        /// <summary>
        /// SnpEffコマンドの引数に変換する。
        /// </summary>
        /// <returns>SnpEff引数</returns>
        internal string ToSnpEffArg()
        {
            return HasFile
                ? $"-c {Path}"
                : "";
        }
    }
}
