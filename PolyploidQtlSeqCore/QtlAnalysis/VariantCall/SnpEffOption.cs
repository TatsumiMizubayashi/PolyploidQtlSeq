namespace PolyploidQtlSeqCore.QtlAnalysis.VariantCall
{
    /// <summary>
    /// SnpEffオプション
    /// </summary>
    internal class SnpEffOption
    {
        private static readonly IReadOnlyDictionary<string, string> _toLongNameDictionary = new Dictionary<string, string>()
        {
            [SnpEffMaxHeap.SHORT_NAME] = SnpEffMaxHeap.LONG_NAME,
            [SnpEffMaxHeap.LONG_NAME] = SnpEffMaxHeap.LONG_NAME,

            [SnpEffConfigFile.SHORT_NAME] = SnpEffConfigFile.LONG_NAME,
            [SnpEffConfigFile.LONG_NAME] = SnpEffConfigFile.LONG_NAME,

            [SnpEffDatabase.SHORT_NAME] = SnpEffDatabase.LONG_NAME,
            [SnpEffDatabase.LONG_NAME] = SnpEffDatabase.LONG_NAME
        };

        /// <summary>
        /// LongName変換項目を辞書に追加する。
        /// </summary>
        /// <param name="dictionary">LongName変換辞書</param>
        public static void AddLongNameKeyValuePair(Dictionary<string, string> dictionary)
        {
            foreach (var keyValuePair in _toLongNameDictionary)
            {
                dictionary.Add(keyValuePair.Key, keyValuePair.Value);
            }
        }


        /// <summary>
        /// SnpEffオプションを作成する。
        /// </summary>
        /// <param name="optionValues">オプションの値</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public SnpEffOption(ISnpEffOption optionValues, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            MaxHeap = new SnpEffMaxHeap(optionValues.SnpEffMaxHeap, parameterDictionary, userOptionDictionary);
            ConfigFile = new SnpEffConfigFile(optionValues.SnpEffConfigFile, parameterDictionary, userOptionDictionary);
            Database = new SnpEffDatabase(optionValues.SnpEffDatabaseName, parameterDictionary, userOptionDictionary);
        }

        /// <summary>
        /// MaxHeapサイズを取得する。
        /// </summary>
        public SnpEffMaxHeap MaxHeap { get; }

        /// <summary>
        /// SnpEff.configファイルを取得する。
        /// </summary>
        public SnpEffConfigFile ConfigFile { get; }

        /// <summary>
        /// SnpEffデータベースを取得する。
        /// </summary>
        public SnpEffDatabase Database { get; }

        /// <summary>
        /// SnpEffが実行可能かどうかを取得する。
        /// </summary>
        public bool CanSneEff => !string.IsNullOrEmpty(Database.Value);

        /// <summary>
        /// パラメータファイルに記載する行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        public string[] ToParameterFileLines()
        {
            return new[]
            {
                MaxHeap.ToParameterFileLine(),
                ConfigFile.ToParameterFileLine(),
                Database.ToParameterFileLine()
            };
        }
    }
}
