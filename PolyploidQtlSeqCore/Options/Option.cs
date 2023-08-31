namespace PolyploidQtlSeqCore.Options
{
    /// <summary>
    /// オプションの抽象クラス
    /// </summary>
    public abstract class Option
    {
        /// <summary>
        /// 該当するLongOptionNameを取得する。
        /// 該当しない場合は空文字を返す。
        /// </summary>
        /// <param name="name">名前</param>
        /// <returns>LongName</returns>
        public string GetLongName(string name)
        {
            if (name == GetShortName()) return GetLongName();
            if (name == GetLongName()) return GetLongName();

            return "";
        }

        /// <summary>
        /// パラメータ辞書から値を設定する。
        /// 該当する情報が辞書に無い場合や、オプションスイッチでユーザー指定されている場合は何もしない。
        /// </summary>
        /// <param name="longNameParameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userSetLongNameDictionary">ユーザー指定オプションLongName辞書</param>
        public void SetValue(IReadOnlyDictionary<string, string> longNameParameterDictionary,
            IReadOnlyDictionary<string, bool> userSetLongNameDictionary)
        {
            if (userSetLongNameDictionary.ContainsKey(GetLongName())) return;
            if (!longNameParameterDictionary.TryGetValue(GetLongName(), out var value)) return;

            SetValue(value);
        }

        /// <summary>
        /// Key-Value文字列を取得する。
        /// </summary>
        /// <returns>Key-Value</returns>
        public string GetKeyValue()
        {
            return $"{GetLongName()}\t{GetStringValue()}";
        }

        /// <summary>
        /// データ検証を行う。
        /// </summary>
        /// <returns>データ検証結果</returns>
        public abstract DataValidationResult Validation();

        /// <summary>
        /// オプションスイッチShortNameを取得する。
        /// </summary>
        /// <returns>ShortName</returns>
        protected abstract string GetShortName();

        /// <summary>
        /// オプションスイッチLongNameを取得する。
        /// </summary>
        /// <returns>LongName</returns>
        protected abstract string GetLongName();

        /// <summary>
        /// 値を設定する。
        /// </summary>
        /// <param name="value">文字列値</param>
        protected abstract void SetValue(string value);

        /// <summary>
        /// 値の文字列値を取得する。
        /// </summary>
        /// <returns>文字列値</returns>
        protected abstract string GetStringValue();
    }
}
