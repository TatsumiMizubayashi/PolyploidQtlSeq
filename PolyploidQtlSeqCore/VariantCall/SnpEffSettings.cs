namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// SnpEff設定
    /// </summary>
    internal class SnpEffSettings
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
        [Obsolete("削除予定")]
        public static void AddLongNameKeyValuePair(Dictionary<string, string> dictionary)
        {
            foreach (var keyValuePair in _toLongNameDictionary)
            {
                dictionary.Add(keyValuePair.Key, keyValuePair.Value);
            }
        }


        /// <summary>
        /// SnpEff設定を作成する。
        /// </summary>
        /// <param name="settingValue">オプションの値</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        [Obsolete("削除予定")]
        public SnpEffSettings(ISnpEffSettingValue settingValue, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            MaxHeap = new SnpEffMaxHeap(settingValue.SnpEffMaxHeap, parameterDictionary, userOptionDictionary);
            ConfigFile = new SnpEffConfigFile(settingValue.SnpEffConfigFile, parameterDictionary, userOptionDictionary);
            Database = new SnpEffDatabase(settingValue.SnpEffDatabaseName, parameterDictionary, userOptionDictionary);
        }

        /// <summary>
        /// SnpEff設定を作成する。
        /// </summary>
        /// <param name="settingValue">オプションの値</param>
        public SnpEffSettings(ISnpEffSettingValue settingValue)
        {
            MaxHeap = new SnpEffMaxHeap(settingValue.SnpEffMaxHeap);
            ConfigFile = new SnpEffConfigFile(settingValue.SnpEffConfigFile);
            Database = new SnpEffDatabase(settingValue.SnpEffDatabaseName);
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
        [Obsolete("削除予定")]
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
