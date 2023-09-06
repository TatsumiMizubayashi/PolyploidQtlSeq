namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// SnpEffデータベース
    /// </summary>
    internal class SnpEffDatabase
    {
        /// <summary>
        /// SnpEffデータベースを作成する。
        /// </summary>
        /// <param name="databaseName">データベース名</param>
        public SnpEffDatabase(string databaseName)
        {
            Value = databaseName;
        }

        /// <summary>
        /// SnpEffデータベース名を取得する。
        /// </summary>
        internal string Value { get; }
    }
}
